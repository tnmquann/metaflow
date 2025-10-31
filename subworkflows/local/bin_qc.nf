#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Perform quality control on genomic bins using BUSCO, CheckM, CheckM2, GUNC, and QUAST
 */

// nf-core modules
include { BUSCO_BUSCO                      } from '../../modules/nf-core/busco/busco/main'
include { CHECKM2_DATABASEDOWNLOAD         } from '../../modules/nf-core/checkm2/databasedownload/main'
include { CHECKM_QA                        } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF                 } from '../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT                  } from '../../modules/nf-core/checkm2/predict/main'
include { CSVTK_CONCAT as CONCAT_BINQC_TSV } from '../../modules/nf-core/csvtk/concat/main'
include { GUNC_DOWNLOADDB                  } from '../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN                         } from '../../modules/nf-core/gunc/run/main'
include { GUNC_MERGECHECKM                 } from '../../modules/nf-core/gunc/mergecheckm/main'
include { UNTAR as BUSCO_UNTAR             } from '../../modules/nf-core/untar/main'
include { UNTAR as CHECKM_UNTAR            } from '../../modules/nf-core/untar/main'
include { GUNZIP                           } from '../../modules/nf-core/gunzip/main'

// local modules
include { QUAST_BINS                       } from '../../modules/local/quast/bins/main'
include { QUAST_BINS_SUMMARY               } from '../../modules/local/quast/binssummary/main'
include { ASSEMBLY_STATS } from '../../modules/local/finalize/assemblybased/assembly_stats'

workflow BIN_QC {
    take:
    ch_bins // channel: [ val(meta), [path(bins)] ] - from BINNING and/or BINNING_REFINEMENT

    main:
    ch_qc_summary = []
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_qc_summaries = Channel.empty()

    /*
    ================================
     * Setup databases
    ================================
     */

    // BUSCO database setup
    if (params.busco_db) {
        ch_busco_db = file(params.busco_db, checkIfExists: true)
        
        // Check if database needs to be extracted
        if (ch_busco_db.toString().endsWith('.tar.gz') || ch_busco_db.toString().endsWith('.tgz')) {
            BUSCO_UNTAR([[id: 'busco_db'], ch_busco_db])
            ch_versions = ch_versions.mix(BUSCO_UNTAR.out.versions)
            ch_busco_db = BUSCO_UNTAR.out.untar.map { it[1] }
        }
    } else {
        ch_busco_db = []
    }

    // CheckM database setup
    if (params.checkm_db) {
        ch_checkm_db = file(params.checkm_db, checkIfExists: true)
    } else if (params.binqc_tool == 'checkm') {
        ch_checkm_db = [[id: 'checkm_db'], file(params.checkm_download_url, checkIfExists: true)]
        CHECKM_UNTAR(ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM_UNTAR.out.versions)
        ch_checkm_db = CHECKM_UNTAR.out.untar.map { it[1] }
    } else {
        ch_checkm_db = []
    }

    // CheckM2 database setup
    if (params.checkm2_db) {
        // CheckM2 expects tuple val(dbmeta), path(db)
        ch_checkm2_db = [[id: 'checkm2_db'], file(params.checkm2_db, checkIfExists: true)]
    } else if (params.binqc_tool == 'checkm2') {
        CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_version)
        ch_versions = ch_versions.mix(CHECKM2_DATABASEDOWNLOAD.out.versions)
        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
    } else {
        ch_checkm2_db = [[id: 'checkm2_db'], []]
    }

    // GUNC database setup
    if (params.gunc_db) {
        ch_gunc_db = file(params.gunc_db, checkIfExists: true)
    } else if (params.run_gunc) {
        GUNC_DOWNLOADDB(params.gunc_db_name)
        ch_versions = ch_versions.mix(GUNC_DOWNLOADDB.out.versions)
        ch_gunc_db = GUNC_DOWNLOADDB.out.db
    } else {
        ch_gunc_db = []
    }

    /*
    ================================
     * Prepare bins for QC
    ================================
     */

    // Convert ArrayBag to proper list and prepare for QUAST
    ch_bins_for_quast = ch_bins
        .map { meta, bins ->
            // Flatten ArrayBag to list if needed
            def bin_list = bins instanceof Collection ? bins.flatten() : [bins]
            [meta, bin_list]
        }

    /*
    ================================
     * Run QUAST & Asselbly statistics on bins
    ================================
     */

    QUAST_BINS(ch_bins_for_quast)
    ch_versions = ch_versions.mix(QUAST_BINS.out.versions)

    QUAST_BINS_SUMMARY(QUAST_BINS.out.quast_bin_summaries.map { it[1] }.collect())
    ch_versions = ch_versions.mix(QUAST_BINS_SUMMARY.out.versions)

    ASSEMBLY_STATS(ch_bins_for_quast)
    ch_versions = ch_versions.mix(ASSEMBLY_STATS.out.versions)
    /*
    ================================
     * Run QC tools based on binqc_tool parameter
    ================================
     */

    if (params.binqc_tool == "busco") {
        // Transpose bins for individual processing
        ch_input_bins_for_qc = ch_bins
            .map { meta, bins ->
                def bin_list = bins instanceof Collection ? bins.flatten() : [bins]
                [meta, bin_list]
            }
            .transpose()

        // Prepare database object depending on type
        if (ch_busco_db && ch_busco_db.extension in ['gz', 'tgz']) {
            ch_busco_db_for_process = BUSCO_UNTAR.out.untar.map { it[1] }
        } else if (ch_busco_db && ch_busco_db.isDirectory()) {
            ch_busco_db_for_process = ch_busco_db
        } else {
            ch_busco_db_for_process = []
        }

        BUSCO_BUSCO(ch_input_bins_for_qc, 'genome', params.busco_db_lineage, ch_busco_db_for_process, [], params.busco_clean)
        ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions)

        ch_qc_summaries = BUSCO_BUSCO.out.batch_summary
            .map { _meta, summary -> [[id: 'busco'], summary] }
            .groupTuple()
        ch_multiqc_files = ch_multiqc_files.mix(
            BUSCO_BUSCO.out.short_summaries_txt.map { it[1] }.flatten()
        )
    } else if (params.binqc_tool == "checkm") {
        // Prepare bins and decompress if needed
        ch_bins_for_checkm = ch_bins
            .map { meta, bins ->
                def bin_list = bins instanceof Collection ? bins.flatten() : [bins]
                [meta, bin_list]
            }
            .transpose()
            .branch { meta, bin ->
                compressed: bin.getName().endsWith('.gz')
                uncompressed: true
            }

        // Decompress gzipped bins
        GUNZIP(ch_bins_for_checkm.compressed)
        ch_versions = ch_versions.mix(GUNZIP.out.versions)

        // Combine decompressed and already uncompressed bins
        ch_bins_uncompressed = GUNZIP.out.gunzip
            .mix(ch_bins_for_checkm.uncompressed)
            .groupTuple()

        // Run CheckM lineage workflow with uncompressed bins
        CHECKM_LINEAGEWF(ch_bins_uncompressed, "fa", ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)

        // Prepare input for CHECKM_QA by joining checkm_output and marker_file
        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, ch_checkm_db ?: []]
            }

        // Run CheckM QA with checkm_db staged as coverage_file
        CHECKM_QA(ch_checkmqa_input, [])
        ch_versions = ch_versions.mix(CHECKM_QA.out.versions)

        // Prepare for concatenation - use a common ID for grouping all samples
        ch_qc_summaries = CHECKM_QA.out.output
            .map { meta, summary -> [[id: 'checkm'], summary] }
            .groupTuple()
        
        ch_multiqc_files = ch_multiqc_files.mix(
            CHECKM_QA.out.output.map { it[1] }.flatten()
        )
    } else if (params.binqc_tool == "checkm2") {
        // Chuẩn bị bins cho CheckM2 - nhóm tất cả bins theo sample
        ch_bins_for_checkm2 = ch_bins
            .map { meta, bins ->
                def bin_list = bins instanceof Collection ? bins.flatten() : [bins]
                [meta, bin_list]
            }

        // Chạy CHECKM2_PREDICT với bins và database
        CHECKM2_PREDICT(ch_bins_for_checkm2, ch_checkm2_db)
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)

        // Gom kết quả cho multiqc và summary
        ch_qc_summaries = CHECKM2_PREDICT.out.checkm2_tsv
            .map { _meta, summary -> [[id: 'checkm2'], summary] }
            .groupTuple()
        ch_multiqc_files = ch_multiqc_files.mix(
            CHECKM2_PREDICT.out.checkm2_tsv.map { meta, tsv -> tsv }
        )
    }
    /*
    ================================
     * Run GUNC for chimera detection (optional)
    ================================
     */

    if (params.run_gunc) {
        ch_input_bins_for_gunc = ch_bins
            .map { meta, bins ->
                def bin_list = bins instanceof Collection ? bins.flatten() : [bins]
                [meta, bin_list]
            }
           // .transpose()

        GUNC_RUN(ch_input_bins_for_gunc, ch_gunc_db)
        ch_versions = ch_versions.mix(GUNC_RUN.out.versions)

        // Make sure to keep directory in sync with modules.conf
        GUNC_RUN.out.maxcss_level_tsv
            .map { _meta, gunc_summary -> gunc_summary }
            .collectFile(
                name: "gunc_summary.tsv",
                keepHeader: true,
                storeDir: "${params.outdir}/Binning/QC",
            )

        if (params.binqc_tool == 'checkm') {
            ch_input_to_mergecheckm = GUNC_RUN.out.maxcss_level_tsv.combine(CHECKM_QA.out.output, by: 0)

            GUNC_MERGECHECKM(ch_input_to_mergecheckm)
            ch_versions.mix(GUNC_MERGECHECKM.out.versions)

            // Make sure to keep directory in sync with modules.conf
            GUNC_MERGECHECKM.out.tsv
                .map { _meta, gunc_checkm_summary -> gunc_checkm_summary }
                .collectFile(
                    name: "gunc_checkm_summary.tsv",
                    keepHeader: true,
                    storeDir: "${params.outdir}/Binning/QC",
                )
        }
    }

    /*
    ================================
     * Concatenate QC summaries
    ================================
     */

    if (ch_qc_summaries) {
        CONCAT_BINQC_TSV(ch_qc_summaries, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CONCAT_BINQC_TSV.out.versions)
        ch_qc_summary = CONCAT_BINQC_TSV.out.csv.map { meta, summary -> summary }
    } else {
        ch_qc_summary = Channel.empty()
    }

    emit:
    qc_summary    = ch_qc_summary
    quast_summary = QUAST_BINS_SUMMARY.out.summary
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}