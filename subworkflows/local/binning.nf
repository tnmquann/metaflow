#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Include nf-core modules
include { METABAT2_METABAT2                                            } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS                         } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                                                      } from '../../modules/nf-core/maxbin2/main'
include { GUNZIP as GUNZIP_BINS                                        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS                                      } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ASSEMBLIES_CONCOCT                          } from '../../modules/nf-core/gunzip/main'

// Include local modules
include { CONVERT_DEPTHS                        } from '../../modules/local/convert_depths'
include { ADJUST_MAXBIN2_EXT                    } from '../../modules/local/adjust_maxbin2_ext'
include { SPLIT_FASTA                           } from '../../modules/local/split_fasta'

// Include nf-core subworkflow
include { FASTA_BINNING_CONCOCT                 } from '../nf-core/fasta_binning_concoct/main'

workflow BINNING {
    take:
    assemblies      // channel: [ val(meta), path(contigs) ] - from MEGAHIT/METASPADES
    bam_bai         // channel: [ val(meta), path(bam), path(bai) ] - from BINNING_BAMABUND

    main:
    versions_ch = Channel.empty()
    ch_metabat2_bins = Channel.empty()
    ch_maxbin2_bins = Channel.empty()
    ch_concoct_bins = Channel.empty()
    ch_all_bins = Channel.empty()

    // Group assemblies and BAMs by key (sample_id + assembler)
    ch_assemblies_grouped = assemblies.map { meta, assembly ->
        ["${meta.id}_${meta.assembler}", meta, assembly]
    }

    ch_bam_bai_grouped = bam_bai.map { meta, bam, bai ->
        ["${meta.id}_${meta.assembler}", meta, bam, bai]
    }

    // Combine and create the bundle format expected by MAG-style BINNING
    // Format: [ val(meta), path(assembly), [path(bams)], [path(bais)] ]
    ch_for_binning = ch_assemblies_grouped
        .combine(ch_bam_bai_grouped, by: 0)
        .map { key, assembly_meta, assembly, bam_meta, bam, bai ->
            // Create arrays for bams and bais - this is critical for CONCOCT
            [assembly_meta, assembly, [bam], [bai]]
        }

    // ==============================================
    // MetaBAT2 Binning
    // ==============================================
    if (!params.skip_metabat2) {
        // Calculate depths using jgi_summarize_bam_contig_depths
        ch_summarizedepth_input = ch_for_binning.map { meta, assembly, bams, bais ->
            [meta, bams, bais]
        }

        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(ch_summarizedepth_input)
        versions_ch = versions_ch.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions)

        ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
            .map { meta, depths ->
                def meta_new = meta.clone()
                meta_new.binner = 'MetaBAT2'
                [meta_new, depths]
            }

        // Combine depths back with assemblies
        ch_metabat2_input = ch_for_binning
            .map { meta, assembly, bams, bais ->
                def meta_new = meta.clone()
                meta_new.binner = 'MetaBAT2'
                [meta_new, assembly, bams, bais]
            }
            .join(ch_metabat_depths, by: 0)
            .map { meta, assembly, bams, bais, depths ->
                [meta, assembly, depths]
            }

        METABAT2_METABAT2(ch_metabat2_input)
        versions_ch = versions_ch.mix(METABAT2_METABAT2.out.versions)

        // Process MetaBAT2 bins
        ch_metabat2_bins_raw = METABAT2_METABAT2.out.fasta
            .transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.bin_id = bin.baseName
                [meta_new, bin]
            }

        // Split unbinned contigs
        ch_metabat2_unbinned = METABAT2_METABAT2.out.unbinned
            .map { meta, unbinned ->
                [meta, unbinned]
            }

        SPLIT_FASTA(ch_metabat2_unbinned)
        versions_ch = versions_ch.mix(SPLIT_FASTA.out.versions)

        ch_metabat2_bins = ch_metabat2_bins_raw.mix(
            SPLIT_FASTA.out.unbinned.transpose()
        )
    }

    // ==============================================
    // MaxBin2 Binning
    // ==============================================
    if (!params.skip_maxbin2) {
        // Convert depths for MaxBin2 format
        CONVERT_DEPTHS(ch_metabat2_input)
        versions_ch = versions_ch.mix(CONVERT_DEPTHS.out.versions)

        ch_maxbin2_input = CONVERT_DEPTHS.out.output
            .map { meta, assembly, reads, depth ->
                def meta_new = meta.clone()
                meta_new.binner = 'MaxBin2'
                [meta_new, assembly, reads, depth]
            }

        MAXBIN2(ch_maxbin2_input)
        versions_ch = versions_ch.mix(MAXBIN2.out.versions)

        // Adjust MaxBin2 file extensions
        ch_maxbin2_bins_to_adjust = MAXBIN2.out.binned_fastas
            .map { meta, bins ->
                [meta, bins]
            }

        ADJUST_MAXBIN2_EXT(ch_maxbin2_bins_to_adjust)

        ch_maxbin2_bins = ADJUST_MAXBIN2_EXT.out.renamed_bins
            .transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.bin_id = bin.baseName
                [meta_new, bin]
            }
    }

    // ==============================================
    // CONCOCT Binning
    // ==============================================
    if (!params.skip_concoct) {
        // CRITICAL: CONCOCT requires:
        // 1. Uncompressed FASTA files
        // 2. Sorted BAM files with proper headers
        // 3. Single BAM per assembly (not arrays)
        ch_assemblies_for_gunzip_concoct = ch_for_binning.map { meta, assembly, bams, bais ->
            [meta, assembly]
        }

        GUNZIP_ASSEMBLIES_CONCOCT(ch_assemblies_for_gunzip_concoct)
        versions_ch = versions_ch.mix(GUNZIP_ASSEMBLIES_CONCOCT.out.versions)

        // Match uncompressed assemblies with their BAMs using composite key
        ch_concoct_input = GUNZIP_ASSEMBLIES_CONCOCT.out.gunzip
            .map { meta, assembly -> 
                def new_meta = meta.clone()
                new_meta.binner = 'CONCOCT'
                ["${meta.id}_${meta.assembler}", new_meta, assembly]
            }
            .combine(
                ch_for_binning.map { meta, assembly, bams, bais ->
                    ["${meta.id}_${meta.assembler}", bams[0], bais[0]]
                },
                by: 0
            )
            .map { key, meta, assembly, bam, bai ->
                [meta, assembly, bam, bai]
            }
            .multiMap { meta, assembly, bam, bai -> 
                bins: [meta, assembly]
                bams: [meta, bam, bai]
            }

        // Run CONCOCT with properly formatted inputs
        FASTA_BINNING_CONCOCT(
            ch_concoct_input.bins,
            ch_concoct_input.bams
        )
        versions_ch = versions_ch.mix(FASTA_BINNING_CONCOCT.out.versions)

        ch_concoct_bins = FASTA_BINNING_CONCOCT.out.bins
            .transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.bin_id = bin.baseName
                [meta_new, bin]
            }
    }

    // Combine all bins from all binners
    ch_all_bins = ch_metabat2_bins
        .mix(ch_maxbin2_bins)
        .mix(ch_concoct_bins)

    emit:
    metabat2_bins   = ch_metabat2_bins      // channel: [ val(meta), path(bin.fa.gz) ]
    maxbin2_bins    = ch_maxbin2_bins       // channel: [ val(meta), path(bin.fa.gz) ]
    concoct_bins    = ch_concoct_bins       // channel: [ val(meta), path(bin.fa.gz) ]
    all_bins        = ch_all_bins           // channel: [ val(meta), path(bin.fa.gz) ]
    versions        = versions_ch
}