#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import modules
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_METABAT2 } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_MAXBIN2  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_CONCOCT  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_COMEBIN  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_SEMIBIN  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_VAMB     } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL                                         } from '../../modules/nf-core/dastool/dastool/main.nf'
include { BINETTE_BINETTE                                         } from '../../modules/local/binette/binette/main.nf'
include { RENAME_PREBINREFINE                                     } from '../../modules/local/rename/prebinrefine/main.nf'
include { RENAME_POSTBINREFINE                                    } from '../../modules/local/rename/postbinrefine/main.nf'
include { CHECKM2_DATABASEDOWNLOAD                             } from '../../modules/nf-core/checkm2/databasedownload/main'

workflow BINNING_REFINEMENT {
    take:
    assemblies      // channel: [ val(meta), path(contigs) ]
    all_bins        // channel: [ val(meta), path(bins) ] - from BINNING.out.all_bins

    main:
    ch_versions = Channel.empty()

    // Remove any refinement metadata and group bins by binner
    ch_bins_grouped = all_bins
        .map { meta, bin ->
            def meta_clean = meta.clone()
            meta_clean.remove('refinement')
            meta_clean.remove('bin_id')
            [meta_clean, bin]
        }
        .groupTuple()
        .map { meta, bins ->
            [meta, bins.flatten()]
        }

    // Prepare bins with renamed format
    RENAME_PREBINREFINE(ch_bins_grouped)

    ch_renamed_bins = RENAME_PREBINREFINE.out.renamed_bins
        .branch {
            metabat2: it[0].binner == 'MetaBAT2'
            maxbin2:  it[0].binner == 'MaxBin2'
            concoct:  it[0].binner == 'CONCOCT'
            comebin:  it[0].binner == 'COMEBin'
            semibin:  it[0].binner == 'SemiBin2'
            vamb:     it[0].binner == 'VAMB'
        }

    // Generate contig2bin tables for each binner
    FASTATOCONTIG2BIN_METABAT2(ch_renamed_bins.metabat2, "fa")
    FASTATOCONTIG2BIN_MAXBIN2(ch_renamed_bins.maxbin2, "fa")
    FASTATOCONTIG2BIN_CONCOCT(ch_renamed_bins.concoct, "fa")
    FASTATOCONTIG2BIN_COMEBIN(ch_renamed_bins.comebin, "fa")
    FASTATOCONTIG2BIN_SEMIBIN(ch_renamed_bins.semibin, "fa")
    FASTATOCONTIG2BIN_VAMB(ch_renamed_bins.vamb, "fa")

    ch_versions = ch_versions.mix(FASTATOCONTIG2BIN_METABAT2.out.versions.first())

    // Collect all contig2bin outputs and group by assembly
    ch_fastatocontig2bin = Channel.empty()
    ch_fastatocontig2bin = ch_fastatocontig2bin
        .mix(FASTATOCONTIG2BIN_METABAT2.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_MAXBIN2.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_CONCOCT.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_COMEBIN.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_SEMIBIN.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_VAMB.out.fastatocontig2bin)
        .map { meta, fastatocontig2bin ->
            def meta_clean = meta.clone()
            meta_clean.remove('binner')
            [meta_clean.id + "_" + meta_clean.assembler, meta_clean, fastatocontig2bin]
        }
        .groupTuple(by: 0)
        .map { key, metas, files ->
            [metas[0], files.flatten()]
        }

    // Branch based on refine_tool parameter
    if (params.refine_tool == 'dastool') {
        // Prepare input for DAS Tool
        ch_input_for_dastool = assemblies
            .map { meta, assembly -> 
                [meta.id + "_" + meta.assembler, meta, assembly] 
            }
            .join(
                ch_fastatocontig2bin.map { meta, files -> 
                    [meta.id + "_" + meta.assembler, files] 
                },
                by: 0
            )
            .map { key, meta, assembly, fastatocontig2bin ->
                def meta_new = meta.clone()
                meta_new.binrefine = 'DASTool'
                [meta_new, assembly, fastatocontig2bin, []]
            }

        // Run DAS Tool for bin refinement
        DASTOOL_DASTOOL(ch_input_for_dastool, [])
        ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions)

        // Process DAS Tool outputs
        ch_refined_output = DASTOOL_DASTOOL.out.bins
            .map { meta, bins ->
                def meta_refined = meta.clone()
                meta_refined.binrefine = 'DASTool'
                [meta_refined, bins]
            }
    } else if (params.refine_tool == 'binette') {
        // CheckM2 database setup (matching BIN_QC pattern)
        if (params.checkm2_db) {
            ch_checkm2_db = Channel.value(file(params.checkm2_db, checkIfExists: true))
        } else {
            CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_version)
            ch_versions = ch_versions.mix(CHECKM2_DATABASEDOWNLOAD.out.versions)
            ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database.map { meta, db -> db }
        }

        // Prepare input for Binette
        ch_input_for_binette = assemblies
            .map { meta, assembly -> 
                def key = "${meta.id}_${meta.assembler}"
                [key, meta, assembly] 
            }
            .join(
                ch_fastatocontig2bin.map { meta, files ->
                    def key = "${meta.id}_${meta.assembler}"
                    [key, files]
                }
            )
            .map { key, meta, assembly, contig2bin_files ->
                def meta_new = meta.clone()
                meta_new.binrefine = 'Binette'
                meta_new.input_mode = 'contig2bin_tables'
                [meta_new, assembly, contig2bin_files, []]
            }

        BINETTE_BINETTE(
            ch_input_for_binette,
            ch_checkm2_db
        )
        ch_versions = ch_versions.mix(BINETTE_BINETTE.out.versions)

        // Process Binette outputs
        ch_refined_output = BINETTE_BINETTE.out.bins
            .map { meta, bins ->
                def meta_refined = meta.clone()
                meta_refined.binrefine = 'Binette'
                [meta_refined, bins]
            }
    }

    // Rename refined bins for both tools
    RENAME_POSTBINREFINE(ch_refined_output)

    emit:
    refined_bins    = RENAME_POSTBINREFINE.out.refined_bins
    refined_unbins  = RENAME_POSTBINREFINE.out.refined_unbins
    qc_reports     = RENAME_POSTBINREFINE.out.qc_reports
    versions        = ch_versions
}