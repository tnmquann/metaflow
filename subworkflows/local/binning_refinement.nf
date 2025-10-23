/*
 * Binning refinement using DAS Tool
 */

include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_METABAT2 } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_MAXBIN2  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as FASTATOCONTIG2BIN_CONCOCT  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL                                          } from '../../modules/nf-core/dastool/dastool/main.nf'
include { RENAME_PREBINREFINE                                      } from '../../modules/local/rename/prebinrefine/main.nf'
include { RENAME_POSTBINREFINE                                     } from '../../modules/local/rename/postbinrefine/main.nf'

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

    // Prepare bins with renamed format for DAS Tool
    RENAME_PREBINREFINE(ch_bins_grouped)

    ch_renamed_bins = RENAME_PREBINREFINE.out.renamed_bins
        .branch {
            metabat2: it[0].binner == 'MetaBAT2'
            maxbin2:  it[0].binner == 'MaxBin2'
            concoct:  it[0].binner == 'CONCOCT'
        }

    // Generate DAS Tool auxiliary files for each binner
    FASTATOCONTIG2BIN_METABAT2(ch_renamed_bins.metabat2, "fa")
    FASTATOCONTIG2BIN_MAXBIN2(ch_renamed_bins.maxbin2, "fa")
    FASTATOCONTIG2BIN_CONCOCT(ch_renamed_bins.concoct, "fa")

    ch_versions = ch_versions.mix(FASTATOCONTIG2BIN_METABAT2.out.versions.first())
    ch_versions = ch_versions.mix(FASTATOCONTIG2BIN_MAXBIN2.out.versions.first())
    ch_versions = ch_versions.mix(FASTATOCONTIG2BIN_CONCOCT.out.versions.first())

    // Collect all fastatocontig2bin outputs and group by assembly
    ch_fastatocontig2bin_for_dastool = Channel.empty()
    ch_fastatocontig2bin_for_dastool = ch_fastatocontig2bin_for_dastool
        .mix(FASTATOCONTIG2BIN_METABAT2.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_MAXBIN2.out.fastatocontig2bin)
        .mix(FASTATOCONTIG2BIN_CONCOCT.out.fastatocontig2bin)
        .map { meta, fastatocontig2bin ->
            def meta_clean = meta.clone()
            meta_clean.remove('binner')
            [meta_clean.id + "_" + meta_clean.assembler, meta_clean, fastatocontig2bin]
        }
        .groupTuple(by: 0)
        .map { key, metas, files ->
            [metas[0], files.flatten()]
        }

    // Prepare input for DAS Tool by joining assemblies with fastatocontig2bin files
    // DASTOOL_DASTOOL signature: input: tuple val(meta), path(contigs), path(bins), path(proteins)
    //                            params: path(db_directory)
    ch_input_for_dastool = assemblies
        .map { meta, assembly -> 
            [meta.id + "_" + meta.assembler, meta, assembly] 
        }
        .join(
            ch_fastatocontig2bin_for_dastool.map { meta, files -> 
                [meta.id + "_" + meta.assembler, files] 
            },
            by: 0
        )
        .map { key, meta, assembly, fastatocontig2bin ->
            // DASTOOL_DASTOOL expects: tuple val(meta), path(contigs), path(bins), path(proteins)
            // Then separate parameter: path(db_directory)
            [meta, assembly, fastatocontig2bin, []]  // Empty list for proteins
        }

    // Run DAS Tool for bin refinement
    // Call signature: DASTOOL_DASTOOL(input_tuple, db_directory)
    DASTOOL_DASTOOL(ch_input_for_dastool, [])
    ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions.first())

    // Process DAS Tool output bins
    ch_dastool_output = DASTOOL_DASTOOL.out.bins
        .map { meta, bins ->
            def meta_refined = meta.clone()
            meta_refined.binner = 'DASTool'
            meta_refined.refinement = 'dastool_refined'
            [meta_refined, bins]
        }

    // Rename refined bins
    RENAME_POSTBINREFINE(ch_dastool_output)

    // Prepare refined bins output (exclude unbinned)
    ch_refined_bins = RENAME_POSTBINREFINE.out.refined_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    // Prepare refined unbinned output
    ch_refined_unbins = RENAME_POSTBINREFINE.out.refined_unbins
        .map { meta, unbinned ->
            def meta_new = meta.clone()
            meta_new.refinement = 'dastool_refined_unbinned'
            meta_new.bin_id = unbinned.baseName
            [meta_new, unbinned]
        }

    emit:
    refined_bins    = ch_refined_bins       // channel: [ val(meta), path(bin) ]
    refined_unbins  = ch_refined_unbins     // channel: [ val(meta), path(unbinned) ]
    versions        = ch_versions           // channel: [ path(versions.yml) ]
}