#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Module imports
include { MERGE_PAIREDENDSEQS } from './modules/local/merge/main'
include { SOURMASH_MANYSKETCH } from './modules/local/sourmash/manysketch/main'
include { SOURMASH_FASTMULTIGATHER } from './modules/local/sourmash/fastmultigather/main'
include { YACHT_RUN } from './modules/local/yacht/run/main'
include { PROCESS_READBASED_RESULTS } from './modules/local/finalize/readbased/process_readbased'
include { CLEANUP } from './modules/local/cleanup/main'
include { RGI_PREPARECARDDB } from './modules/local/rgi/preparecarddb/main'
include { RGI_BWT } from './modules/local/rgi/bwt/main'

// Subworkflow imports
include { PREPROCESS } from './subworkflows/local/preprocess'
include { UTILS_NFSCHEMA_PLUGIN } from './subworkflows/nf-core/utils_nfschema_plugin/main'

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Function to create input channel from a CSV file
def createCsvInputChannel(input_path) {
    return Channel
        .fromPath(input_path)
        .splitCsv(header:true, sep:',', strip:true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id ?: "sample_${System.currentTimeMillis()}"
            meta.run_id = row.run_id ?: 'default_run'
            meta.group = row.group ?: 'default_group'
            meta.single_end = false // Assuming paired-end as per your setup

            // Validate read files exist
            def reads = []
            if (row.short_reads_1?.trim() && row.short_reads_2?.trim()) {
                def read1_path = row.short_reads_1.trim()
                def read2_path = row.short_reads_2.trim()
                if (!file(read1_path).exists()) {
                    exit 1, "Read file 1 does not exist for sample ${meta.id}: ${read1_path}"
                }
                if (!file(read2_path).exists()) {
                    exit 1, "Read file 2 does not exist for sample ${meta.id}: ${read2_path}"
                }
                reads = [file(read1_path), file(read2_path)]
            } else {
                exit 1, "Missing or invalid paired-end read files for sample: ${meta.id}. Both short_reads_1 and short_reads_2 must be provided."
            }
            return tuple(meta, reads)
        }
}

workflow {
    // Validate parameters and log summary
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Create input channel based on the specified format
    if (params.input_format == 'csv') {
        input_ch = createCsvInputChannel(params.input)
    } else if (params.input_format == 'directory') {
        // Use fromFilePairs for directory input to find paired-end FASTQ files
        input_ch = Channel.fromFilePairs("${params.input}/*_{1,2}*.fastq.gz")
            .map { sample_id, reads ->
                def meta = [:]
                meta.id = sample_id
                meta.single_end = false
                // Note: run_id and group are not available from filenames, using defaults
                meta.run_id = 'default_run'
                meta.group = 'default_group'
                return tuple(meta, reads)
            }
    } else {
        exit 1, "Invalid input_format: '${params.input_format}'. Must be 'csv' or 'directory'."
    }

    // Run preprocessing subworkflow including QC, trimming, and host removal
    PREPROCESS(input_ch)

    // Create a fork in the workflow after PREPROCESS by duplicating the channel
    def cleaned_reads_ch = PREPROCESS.out.cleaned_reads

    // RGI Branch - Only execute if enable_rgi_bwt is true
    if (params.enable_rgi_bwt) {
        // Transform reads for RGI_BWT
        def rgi_reads = cleaned_reads_ch
            .map { meta, reads -> 
                // Split the reads array into individual paths for RGI
                def read1 = reads[0]
                def read2 = reads[1]
                return [meta, read1, read2]
            }

        // Prepare RGI database if not provided
        ch_rgi_db = params.rgi_preparecarddb_dir ? 
            Channel.value(file(params.rgi_preparecarddb_dir)) : 
            RGI_PREPARECARDDB([]).db

        // Run RGI BWT using cleaned reads and prepared database
        RGI_BWT(
            rgi_reads,
            ch_rgi_db
        )
    }

    // MERGE Branch - existing workflow continues
    def merge_reads = cleaned_reads_ch
    MERGE_PAIREDENDSEQS(merge_reads)

    // Prepare input for SOURMASH_MANYSKETCH
    ch_manysketch_csv_input = MERGE_PAIREDENDSEQS.out.merged_seqs
        .collectFile(
            name: "${params.outdir}/Sourmash - YACHT/manysketch_manifest.csv",
            newLine: true,
            seed: "name,genome_filename,protein_filename",  // Header line
            storeDir: "${params.outdir}/Sourmash - YACHT"
        ) { meta, merged_file ->
            // Each line will be: sampleId,merged_fastq_path,
            "${meta.id},${merged_file},"
        }
        .map { csv_file -> 
            // Return the CSV file for SOURMASH_MANYSKETCH
            csv_file
        }


    // Run taxonomic profiling with Sourmash
    SOURMASH_MANYSKETCH(
        ch_manysketch_csv_input
    )

    // SOURMASH_FASTMULTIGATHER expects: tuple val(meta), path(manysketch_zip)
    // SOURMASH_MANYSKETCH emits: path("batch.manysketch.zip")
    // We need to re-introduce a meta map if subsequent processes need it.
    // For now, let's assume a single batch operation for fastmultigather.
    // If per-sample fastmultigather is needed, SOURMASH_MANYSKETCH would need to emit per-sample zips.
    ch_sketch_zip_with_meta = SOURMASH_MANYSKETCH.out.sketch_zip_file.map{ file -> [[id:'batch'], file] }

    SOURMASH_FASTMULTIGATHER(
        ch_sketch_zip_with_meta,
        file(params.sourmash_database, checkIfExists: true)
    )

    // YACHT_RUN expects: tuple val(meta), path(zip_files_dir)
    // SOURMASH_MANYSKETCH emits: path("manysketch_output/zip_files"), which is a directory
    ch_zip_files_dir_with_meta = SOURMASH_MANYSKETCH.out.zip_files_dir.map{ dir -> [[id:'batch'], dir] }

    YACHT_RUN(
        ch_zip_files_dir_with_meta,
        file(params.yacht_database, checkIfExists: true)
    )

    // Process and combine results
    // PROCESS_READBASED_RESULTS expects:
    // tuple val(meta), path(gather_csv)
    // tuple val(meta), path(yacht_xlsx)
    // SOURMASH_FASTMULTIGATHER emits: tuple val(meta), path("fastmultigather/*_sourmash_gather.csv")
    // YACHT_RUN emits: tuple val(meta), path("yacht_results/*.xlsx")
    // These should align if the meta map has 'id':'batch'

    PROCESS_READBASED_RESULTS(
        SOURMASH_FASTMULTIGATHER.out.gather_csv,
        YACHT_RUN.out.yacht_xlsx
    )

    // Optional cleanup
    if (params.cleanup) {
        CLEANUP(
            PROCESS_READBASED_RESULTS.out.final_results.collect() // collect might be needed if it's a channel
        )
    }
}