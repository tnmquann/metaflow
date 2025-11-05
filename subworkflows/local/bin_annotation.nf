#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Annotate genomic bins using Bakta or Prokka
 */

// Import nf-core modules for bin annotation
include { BAKTA_BAKTADBDOWNLOAD          } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA as BAKTA_BINS      } from '../../modules/nf-core/bakta/bakta/main'
include { PROKKA as PROKKA_BINS          } from '../../modules/nf-core/prokka/main'

workflow BIN_ANNOTATION {
    take:
    ch_bins // channel: [ val(meta), path(bin) ] - from BINNING and/or BINNING_REFINEMENT

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Prepare bins for annotation - ensure single file per tuple
    ch_bins_for_annotation = ch_bins
        .map { meta, bin ->
            // Ensure bin is a single file, not a list
            def bin_file = bin instanceof List ? bin[0] : bin
            // Create new meta with bin-specific ID
            def meta_new = meta.clone()
            meta_new.id = bin_file.baseName
            meta_new.bin_id = bin_file.baseName
            [meta_new, bin_file]
        }

    /*
    ================================
     * Setup annotation databases
    ================================
     */

    if (params.bin_annotation_tool == 'bakta') {
        // BAKTA database setup
        if (params.annotation_bakta_db) {
            ch_bakta_db = Channel
                .fromPath(params.annotation_bakta_db, checkIfExists: true)
                .first()
        } else {
            BAKTA_BAKTADBDOWNLOAD(params.annotation_bakta_db_downloadtype)
            ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
            ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        }

        // Run BAKTA annotation on bins
        BAKTA_BINS(
            ch_bins_for_annotation,
            ch_bakta_db,
            [],  // proteins
            [],  // prodigal_tf
            [],  // regions
            []   // hmms
        )
        ch_versions = ch_versions.mix(BAKTA_BINS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            BAKTA_BINS.out.txt.map { it[1] }.flatten()
        )

        // Output channels
        ch_annotation_faa = BAKTA_BINS.out.faa
        ch_annotation_fna = BAKTA_BINS.out.fna
        ch_annotation_gbk = BAKTA_BINS.out.gbff
        ch_annotation_gff = BAKTA_BINS.out.gff
        ch_annotation_tsv = BAKTA_BINS.out.tsv

    } else if (params.bin_annotation_tool == 'prokka') {
        // Run PROKKA annotation on bins
        PROKKA_BINS(
            ch_bins_for_annotation,
            [],  // proteins
            []   // prodigal_tf
        )
        ch_versions = ch_versions.mix(PROKKA_BINS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            PROKKA_BINS.out.txt.map { it[1] }.flatten()
        )

        // Output channels
        ch_annotation_faa = PROKKA_BINS.out.faa
        ch_annotation_fna = PROKKA_BINS.out.fna
        ch_annotation_gbk = PROKKA_BINS.out.gbk
        ch_annotation_gff = PROKKA_BINS.out.gff
        ch_annotation_tsv = PROKKA_BINS.out.tsv
    }

    emit:
    faa           = ch_annotation_faa        // channel: [ val(meta), path(faa) ]
    fna           = ch_annotation_fna        // channel: [ val(meta), path(fna) ]
    gbk           = ch_annotation_gbk        // channel: [ val(meta), path(gbk) ]
    gff           = ch_annotation_gff        // channel: [ val(meta), path(gff) ]
    tsv           = ch_annotation_tsv        // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files         // channel: [ path(multiqc_files) ]
    versions      = ch_versions              // channel: [ path(versions.yml) ]
}

