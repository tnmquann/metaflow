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
    assemblies      // channel: [ val(meta), path(contigs) ]
    bam_bai         // channel: [ val(meta), [path(bams)], [path(bais)] ] - from BINNING_BAMABUND

    main:
    versions_ch = Channel.empty()
    ch_metabat2_bins = Channel.empty()
    ch_maxbin2_bins = Channel.empty()
    ch_concoct_bins = Channel.empty()
    ch_all_bins = Channel.empty()

    // Combine assemblies with their BAM/BAI files
    // Format: [ val(meta), path(assembly), [path(bams)], [path(bais)] ]
    ch_for_binning = assemblies
        .map { meta, assembly -> [meta.id + "_" + meta.assembler, meta, assembly] }
        .join(
            bam_bai.map { meta, bams, bais -> [meta.id + "_" + meta.assembler, bams, bais] },
            by: 0
        )
        .map { key, meta, assembly, bams, bais ->
            // Ensure assembly is single file, bams/bais are arrays
            [meta, assembly, bams, bais]
        }

    // ==============================================
    // MetaBAT2 Binning
    // ==============================================
    if (!params.skip_metabat2) {
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

        ch_maxbin2_bins_to_adjust = MAXBIN2.out.binned_fastas

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
    // CONCOCT Binning - Fixed to match MAG pattern
    // ==============================================
    if (!params.skip_concoct) {
        // Step 1: Decompress assemblies for CONCOCT (requires uncompressed FASTA)
        ch_assemblies_for_gunzip = ch_for_binning.map { meta, assembly, bams, bais ->
            def meta_new = meta.clone()
            meta_new.binner = 'CONCOCT'
            [meta_new, assembly]
        }

        GUNZIP_ASSEMBLIES_CONCOCT(ch_assemblies_for_gunzip)
        versions_ch = versions_ch.mix(GUNZIP_ASSEMBLIES_CONCOCT.out.versions)

        // Step 2: Prepare multiMap input matching the FASTA_BINNING_CONCOCT signature
        // CONCOCT expects: bins channel [meta, fasta] and bams channel [meta, [bams], [bais]]
        ch_concoct_input = ch_for_binning
            .map { meta, assembly, bams, bais ->
                def meta_new = meta.clone()
                meta_new.binner = 'CONCOCT'
                [meta.id + "_" + meta.assembler, meta_new, assembly, bams, bais]
            }
            .join(
                GUNZIP_ASSEMBLIES_CONCOCT.out.gunzip.map { meta, fasta ->
                    [meta.id + "_" + meta.assembler, fasta]
                },
                by: 0
            )
            .map { key, meta, assembly_gz, bams, bais, fasta_uncompressed ->
                [meta, fasta_uncompressed, bams, bais]
            }
            .multiMap { meta, fasta, bams, bais ->
                bins: [meta, fasta]
                bams: [meta, bams, bais]
            }

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
    metabat2_bins   = ch_metabat2_bins
    maxbin2_bins    = ch_maxbin2_bins
    concoct_bins    = ch_concoct_bins
    all_bins        = ch_all_bins
    versions        = versions_ch
}