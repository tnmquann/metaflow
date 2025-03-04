#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERGE_SEQUENCES } from './modules/merge_sequences'
include { SOURMASH_MANYSKETCH } from './modules/sourmash/manysketch'
include { SOURMASH_FASTMULTIGATHER } from './modules/sourmash/fastmultigather'
include { YACHT } from './modules/yacht/yacht'
include { PROCESS_RESULTS } from './modules/process_results'
// include { CLEANUP } from './modules/cleanup'

workflow {
    merged_seqs = MERGE_SEQUENCES(params.trimmed_fastq)
    manysketch_results = SOURMASH_MANYSKETCH(merged_seqs)
    sourmash_results = SOURMASH_FASTMULTIGATHER(
        manysketch_results.sketch_zip,
        params.sourmash_database,
        merged_seqs
    )
    yacht_results = YACHT(manysketch_results.zip_files, params.yacht_database)
    
    PROCESS_RESULTS(
        sourmash_results,     // this will be the text_files output
        yacht_results,
        merged_seqs,
        manysketch_results.manysketch_dir,
        manysketch_results.sketch_zip
    )
    
    // // Add cleanup as the final step
    // if (params.cleanup) {
    //     CLEANUP(PROCESS_RESULTS.out)
    // }
}

workflow SOURMASH_ONLY {
    merged = MERGE_SEQUENCES(params.trimmed_fastq)
    def manyRes = MANYSKETCH(merged.out)
    FASTMULTIGATHER(manyRes.out.sketch_zip, params.sourmash_database, merged.out)
}

workflow YACHT_ONLY {
    MANYSKETCH_RESULT = SOURMASH_ONLY()
    YACHT(MANYSKETCH_RESULT.out.zip_files, params.yacht_database)
}