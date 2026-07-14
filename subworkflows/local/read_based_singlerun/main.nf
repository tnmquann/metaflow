#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules
include { SOURMASH_SKETCH as SOURMASH_SKETCH_META } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER as SOURMASH_GATHER_META } from '../../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_TAXANNOTATE as SOURMASH_TAXANNOTATE_META } from '../../../modules/nf-core/sourmash/taxannotate/main'

// Import local modules
include { MERGE_PAIREDENDSEQS } from '../../../modules/local/merge/main'
include { SOURMASH_SIG_TO_ZIP } from '../../../modules/local/sourmash/sig_to_zip/main'
include { SOURMASH_TAXMETAGENOME as SOURMASH_TAXMETAGENOME_SINGLESKETCH } from '../../../modules/local/sourmash/taxmetagenome/main'
include { YACHT_RUN_SINGLESKETCH } from '../../../modules/local/yacht/run_singlesketch/main'
include { PROCESS_READBASED_RESULTS_SINGLESKETCH } from '../../../modules/local/finalize/process_readbased_singlesketch/main'
include { POSTPROCESS_READBASED } from '../../../modules/local/finalize/postprocess_readbased/main'
include { MERGE_READBASED_POSTPROCESS_INPUTS } from '../../../modules/local/finalize/merge_readbased_postprocess_inputs/main'
include { RGI_PREPARECARDDB } from '../../../modules/local/rgi/preparecarddb/main'
include { RGI_BWT } from '../../../modules/local/rgi/bwt/main'

def deriveReadbasedBatchPrefix(outdir) {
    def normalized = outdir?.toString()?.trim()?.replaceAll('/+$', '')
    if (!normalized) {
        throw new IllegalArgumentException('Cannot derive a read-based post-processing prefix from an empty --outdir value.')
    }

    def basename = normalized.tokenize('/').last()
    def prefix = basename
        .replaceAll('[^A-Za-z0-9._-]+', '_')
        .replaceAll('^[._-]+', '')
        .replaceAll('[._-]+$', '')

    if (!prefix) {
        throw new IllegalArgumentException("Cannot derive a filesystem-safe read-based post-processing prefix from --outdir '${outdir}'.")
    }
    return prefix
}

def deriveReadbasedSinglePrefix(meta, input_format) {
    def raw_prefix = meta?.id?.toString()?.trim()
    if (!raw_prefix) {
        throw new IllegalArgumentException('Cannot derive a read-based post-processing prefix because the sample metadata has no id.')
    }

    if (input_format == 'directory') {
        raw_prefix = raw_prefix.replaceFirst('([._-])R?[12]$', '')
    }

    def prefix = raw_prefix
        .replaceAll('[^A-Za-z0-9._-]+', '_')
        .replaceAll('^[._-]+', '')
        .replaceAll('[._-]+$', '')

    if (!prefix) {
        throw new IllegalArgumentException("Cannot derive a filesystem-safe read-based post-processing prefix from sample id '${meta.id}'.")
    }
    return prefix
}

workflow READ_BASED_SINGLERUN {
    take:
    cleaned_reads_ch // Channel: [ val(meta), [ reads ] ]

    main:
    versions_ch = channel.empty()
    ch_rgi_results = channel.empty() // Channel for RGI results if enabled
    ch_postprocess_default = channel.empty()
    ch_postprocess_phyloseq = channel.empty()
    ch_postprocess_taxburst = channel.empty()
    ch_postprocess_rgi_bwt = channel.empty()

    // Validate required parameters for this subworkflow
    if (!params.sourmash_database) {
        error "Missing required parameter for read-based single-run profiling: --sourmash_database"
    }
    if (!params.sourmash_taxonomy_csv) {
        error "Missing required parameter for read-based single-run profiling: --sourmash_taxonomy_csv (required by SOURMASH_TAXANNOTATE)"
    }

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
            ch_rgi_db = channel.value(file(params.rgi_preparecarddb_dir))
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
    versions_ch = versions_ch.mix(SOURMASH_SKETCH_META.out.versions_sourmash.first())

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
    versions_ch = versions_ch.mix(SOURMASH_GATHER_META.out.versions_sourmash.first())

    // Step 4: Run sourmash tax annotate (per sample)
    SOURMASH_TAXANNOTATE_META(
        SOURMASH_GATHER_META.out.result,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_TAXANNOTATE_META.out.versions_sourmash.first())

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

        process_singlesketch_script = file("${projectDir}/bin/py_scripts/read_based/process_singlesketch_results.py", checkIfExists: true)

        PROCESS_READBASED_RESULTS_SINGLESKETCH(
            ch_combined_for_processing.map { meta, gather_csv, _yacht_xlsx -> [meta, gather_csv] },
            ch_combined_for_processing.map { meta, _gather_csv, yacht_xlsx -> [meta, yacht_xlsx] },
            process_singlesketch_script
        )
        versions_ch = versions_ch.mix(PROCESS_READBASED_RESULTS_SINGLESKETCH.out.versions.first())

        if (params.readbased_postprocess) {
            def postprocess_script = file("${projectDir}/bin/py_scripts/read_based/post_processing.py", checkIfExists: true)
            def postprocess_base_ch = PROCESS_READBASED_RESULTS_SINGLESKETCH.out.final_results
                .map { meta, final_results_dir ->
                    def output_prefix = deriveReadbasedSinglePrefix(meta, params.input_format)
                    def merged_sourmash_yacht = file("${final_results_dir}/${meta.id}_merged_sourmash_yacht.csv", checkIfExists: true)
                    [
                        meta,
                        output_prefix,
                        params.readbased_postprocess,
                        true,
                        true,
                        params.postprocess_options ?: ' ',
                        merged_sourmash_yacht
                    ]
                }

            def postprocess_input_ch
            if (params.enable_rgi_bwt) {
                postprocess_input_ch = postprocess_base_ch
                    .join(RGI_BWT.out.outdir, by: 0)
                    .map { meta, output_prefix, mode, single_sample_layout, rgi_gene_mapping_only, postprocess_options, merged_sourmash_yacht, rgi_dir ->
                        [meta, output_prefix, mode, single_sample_layout, rgi_gene_mapping_only, postprocess_options, merged_sourmash_yacht, [rgi_dir]]
                    }
            } else {
                postprocess_input_ch = postprocess_base_ch
                    .map { meta, output_prefix, mode, single_sample_layout, rgi_gene_mapping_only, postprocess_options, merged_sourmash_yacht ->
                        [meta, output_prefix, mode, single_sample_layout, rgi_gene_mapping_only, postprocess_options, merged_sourmash_yacht, []]
                    }
            }

            def aggregate_output_prefix = deriveReadbasedBatchPrefix(params.outdir)
            def aggregate_meta = [id: aggregate_output_prefix]
            def aggregate_merged_inputs_ch = PROCESS_READBASED_RESULTS_SINGLESKETCH.out.final_results
                .map { meta, final_results_dir ->
                    file("${final_results_dir}/${meta.id}_merged_sourmash_yacht.csv", checkIfExists: true)
                }
                .collect()
                .map { merged_inputs -> [aggregate_meta, merged_inputs] }

            MERGE_READBASED_POSTPROCESS_INPUTS(aggregate_merged_inputs_ch)
            versions_ch = versions_ch.mix(MERGE_READBASED_POSTPROCESS_INPUTS.out.versions)

            def aggregate_rgi_bundle_ch = params.enable_rgi_bwt ?
                RGI_BWT.out.outdir
                    .map { _meta, rgi_dir -> rgi_dir }
                    .collect()
                    .map { rgi_dirs -> [rgi_dirs] } :
                channel.value([[]])

            def aggregate_postprocess_input_ch = MERGE_READBASED_POSTPROCESS_INPUTS.out.merged_csv
                .combine(aggregate_rgi_bundle_ch)
                .map { meta, merged_sourmash_yacht, rgi_dirs ->
                    [
                        meta,
                        aggregate_output_prefix,
                        params.readbased_postprocess,
                        false,
                        true,
                        params.postprocess_options ?: ' ',
                        merged_sourmash_yacht,
                        rgi_dirs
                    ]
                }

            POSTPROCESS_READBASED(postprocess_input_ch.mix(aggregate_postprocess_input_ch), postprocess_script)
            versions_ch = versions_ch.mix(POSTPROCESS_READBASED.out.versions)
            ch_postprocess_default = POSTPROCESS_READBASED.out.default_output
            ch_postprocess_phyloseq = POSTPROCESS_READBASED.out.phyloseq
            ch_postprocess_taxburst = POSTPROCESS_READBASED.out.taxburst
            ch_postprocess_rgi_bwt = POSTPROCESS_READBASED.out.rgi_bwt
        }
    }

    emit:
    versions = versions_ch.ifEmpty(null)
    results  = params.skip_yacht ? channel.empty() : PROCESS_READBASED_RESULTS_SINGLESKETCH.out.final_results

    rgi_results = ch_rgi_results
    postprocess_default = ch_postprocess_default
    postprocess_phyloseq = ch_postprocess_phyloseq
    postprocess_taxburst = ch_postprocess_taxburst
    postprocess_rgi_bwt = ch_postprocess_rgi_bwt
    gather_csv  = SOURMASH_GATHER_META.out.result
    taxannotate = SOURMASH_TAXANNOTATE_META.out.result
    metagenome_classification = SOURMASH_TAXMETAGENOME_SINGLESKETCH.out.genome_classification
    single_sketches = SOURMASH_SKETCH_META.out.signatures
    single_sketches_zip = SOURMASH_SIG_TO_ZIP.out.sig_zip

    yacht_results = params.skip_yacht ? channel.empty() : YACHT_RUN_SINGLESKETCH.out.yacht_xlsx
}
