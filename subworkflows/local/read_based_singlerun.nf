#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules
include { SOURMASH_SKETCH as SOURMASH_SKETCH_META } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER as SOURMASH_GATHER_META } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_TAXANNOTATE as SOURMASH_TAXANNOTATE_META } from '../../modules/nf-core/sourmash/taxannotate/main'

// Import local modules
include { MERGE_PAIREDENDSEQS } from '../../modules/local/merge/main'
include { SOURMASH_SIG_TO_ZIP } from '../../modules/local/sourmash/sig_to_zip/main'
include { SOURMASH_TAXMETAGENOME as SOURMASH_TAXMETAGENOME_SINGLESKETCH } from '../../modules/local/sourmash/taxmetagenome/main'
include { YACHT_RUN_SINGLESKETCH } from '../../modules/local/yacht/run_singlesketch/main'
include { PROCESS_READBASED_RESULTS_SINGLESKETCH } from '../../modules/local/finalize/readbased/process_readbased_singlesketch'
include { RGI_PREPARECARDDB } from '../../modules/local/rgi/preparecarddb/main'
include { RGI_BWT } from '../../modules/local/rgi/bwt/main'

workflow READ_BASED_SINGLERUN {
    take:
    cleaned_reads_ch // Channel: [ val(meta), [ reads ] ]

    main:
    versions_ch = Channel.empty()
    ch_rgi_results = Channel.empty() // Channel for RGI results if enabled

    // ============================================================
    // RGI Branch - Keep existing logic unchanged
    // ============================================================
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

    // ============================================================
    // SOURMASH/YACHT Branch - Single sample processing
    // Each sample is processed individually through the pipeline
    // ============================================================
    
    // Step 1: Merge paired-end reads (per sample)
    MERGE_PAIREDENDSEQS(cleaned_reads_ch)
    versions_ch = versions_ch.mix(MERGE_PAIREDENDSEQS.out.versions)

    // Step 2: Run sourmash sketch (per sample)
    // Input: merged sequence from each sample
    SOURMASH_SKETCH_META(MERGE_PAIREDENDSEQS.out.merged_seqs)
    versions_ch = versions_ch.mix(SOURMASH_SKETCH_META.out.versions.first())

    // Step 2b: Convert .sig to .sig.zip for YACHT compatibility
    SOURMASH_SIG_TO_ZIP(SOURMASH_SKETCH_META.out.signatures)
    versions_ch = versions_ch.mix(SOURMASH_SIG_TO_ZIP.out.versions.first())

    // Step 3: Run sourmash gather (per sample)
    // Input: signature from each sample + database
    SOURMASH_GATHER_META(
        SOURMASH_SKETCH_META.out.signatures,
        file(params.sourmash_database, checkIfExists: true),
        false,  // save_unassigned
        false,  // save_matches_sig
        false,  // save_prefetch
        false   // save_prefetch_csv
    )
    versions_ch = versions_ch.mix(SOURMASH_GATHER_META.out.versions.first())

    // Step 4: Run sourmash tax annotate (per sample)
    SOURMASH_TAXANNOTATE_META(
        SOURMASH_GATHER_META.out.result,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_TAXANNOTATE_META.out.versions.first())

    // Step 5: Run sourmash tax metagenome (per sample)
    SOURMASH_TAXMETAGENOME_SINGLESKETCH(
        SOURMASH_TAXANNOTATE_META.out.result,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_TAXMETAGENOME_SINGLESKETCH.out.versions.first())

    // Step 6: YACHT branch (parallel with step 4-5)
    // Only run YACHT + PROCESS_READBASED_RESULTS_SINGLESKETCH when not skipping
    if (!params.skip_yacht) {
        // Use .sig.zip files for YACHT
        YACHT_RUN_SINGLESKETCH(
            SOURMASH_SIG_TO_ZIP.out.sig_zip,
            file(params.yacht_database, checkIfExists: true)
        )
        versions_ch = versions_ch.mix(YACHT_RUN_SINGLESKETCH.out.versions.first())

        // Join tax annotate results with YACHT results for processing
        ch_taxannotate_for_processing = SOURMASH_TAXANNOTATE_META.out.result
        ch_yacht_results = YACHT_RUN_SINGLESKETCH.out.yacht_xlsx

        // Join channels by meta to process together
        ch_combined_for_processing = ch_taxannotate_for_processing
            .join(ch_yacht_results, by: 0)
            .map { meta, gather_csv, yacht_xlsx ->
                [meta, gather_csv, yacht_xlsx]
            }

        PROCESS_READBASED_RESULTS_SINGLESKETCH(
            ch_combined_for_processing.map { meta, gather_csv, yacht_xlsx -> [meta, gather_csv] },
            ch_combined_for_processing.map { meta, gather_csv, yacht_xlsx -> [meta, yacht_xlsx] }
        )
        versions_ch = versions_ch.mix(PROCESS_READBASED_RESULTS_SINGLESKETCH.out.versions.first())
    }

    emit:
    versions = versions_ch.ifEmpty(null)
    results  = params.skip_yacht ? Channel.empty() : PROCESS_READBASED_RESULTS_SINGLESKETCH.out.final_results

    rgi_results = ch_rgi_results
    gather_csv  = SOURMASH_GATHER_META.out.result
    taxannotate = SOURMASH_TAXANNOTATE_META.out.result
    metagenome_classification = SOURMASH_TAXMETAGENOME_SINGLESKETCH.out.genome_classification
    single_sketches = SOURMASH_SKETCH_META.out.signatures
    single_sketches_zip = SOURMASH_SIG_TO_ZIP.out.sig_zip

    yacht_results = params.skip_yacht ? Channel.empty() : YACHT_RUN_SINGLESKETCH.out.yacht_xlsx
}