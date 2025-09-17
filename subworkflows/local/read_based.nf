#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import modules
include { MERGE_PAIREDENDSEQS } from '../../modules/local/merge/main'
include { SOURMASH_MANYSKETCH } from '../../modules/local/sourmash/manysketch/main'
include { SOURMASH_FASTMULTIGATHER } from '../../modules/local/sourmash/fastmultigather/main'
include { YACHT_RUN } from '../../modules/local/yacht/run/main'
include { PROCESS_READBASED_RESULTS } from '../../modules/local/finalize/readbased/process_readbased'
include { RGI_PREPARECARDDB } from '../../modules/local/rgi/preparecarddb/main'
include { RGI_BWT } from '../../modules/local/rgi/bwt/main'

workflow READ_BASED {
    take:
    cleaned_reads_ch // Channel: [ val(meta), [ reads ] ]

    main:
    versions_ch = Channel.empty()
    ch_rgi_results = Channel.empty() // Channel for RGI results if enabled

    // RGI Branch - Only execute if enable_rgi_bwt is true 
    if (params.enable_rgi_bwt) {
        // Transform reads for RGI_BWT
        def rgi_reads = cleaned_reads_ch
            .map { meta, reads -> 
                def read1 = reads[0]
                def read2 = reads[1]
                return [meta, read1, read2]
            }

        // Prepare RGI database if not provided
        if (params.rgi_preparecarddb_dir) {
            ch_rgi_db = Channel.value(file(params.rgi_preparecarddb_dir))
        } else {
            RGI_PREPARECARDDB([])
            ch_rgi_db = RGI_PREPARECARDDB.out.db
            versions_ch = versions_ch.mix(RGI_PREPARECARDDB.out.versions)
        }

        // Run RGI BWT
        RGI_BWT(
            rgi_reads,
            ch_rgi_db
        )
        versions_ch = versions_ch.mix(RGI_BWT.out.versions)
        ch_rgi_results = RGI_BWT.out.outdir
    }

    // MERGE Branch
    MERGE_PAIREDENDSEQS(cleaned_reads_ch)
    versions_ch = versions_ch.mix(MERGE_PAIREDENDSEQS.out.versions)

    // Prepare input for SOURMASH_MANYSKETCH
    ch_manysketch_csv_input = MERGE_PAIREDENDSEQS.out.merged_seqs
        .collectFile(
            name: "${params.outdir}/Sourmash - YACHT/manysketch_manifest.csv",
            newLine: true,
            seed: "name,genome_filename,protein_filename",
            storeDir: "${params.outdir}/Sourmash - YACHT"
        ) { meta, merged_file ->
            "${meta.id},${merged_file},"
        }
        .map { csv_file -> csv_file }

    // Run Sourmash pipeline
    SOURMASH_MANYSKETCH(ch_manysketch_csv_input)
    versions_ch = versions_ch.mix(SOURMASH_MANYSKETCH.out.versions)

    ch_sketch_zip_with_meta = SOURMASH_MANYSKETCH.out.sketch_zip_file
        .map{ file -> [[id:'batch'], file] }

    SOURMASH_FASTMULTIGATHER(
        ch_sketch_zip_with_meta,
        file(params.sourmash_database, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_FASTMULTIGATHER.out.versions)

    ch_zip_files_dir_with_meta = SOURMASH_MANYSKETCH.out.zip_files_dir
        .map{ dir -> [[id:'batch'], dir] }

    YACHT_RUN(
        ch_zip_files_dir_with_meta,
        file(params.yacht_database, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(YACHT_RUN.out.versions)

    // Process results
    PROCESS_READBASED_RESULTS(
        SOURMASH_FASTMULTIGATHER.out.gather_csv,
        YACHT_RUN.out.yacht_xlsx
    )
    versions_ch = versions_ch.mix(PROCESS_READBASED_RESULTS.out.versions)

    emit:
    versions = versions_ch.ifEmpty(null)
    results = PROCESS_READBASED_RESULTS.out.final_results
    rgi_results = ch_rgi_results
}