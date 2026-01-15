#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESS } from '../subworkflows/local/preprocess'
include { READ_BASED } from '../subworkflows/local/read_based'
include { READ_BASED_SINGLERUN } from '../subworkflows/local/read_based_singlerun'
include { ASSEMBLY_BASED } from '../subworkflows/local/assembly_based'
include { BINNING_BAMABUND } from '../subworkflows/local/binning_bamabund'
include { CLEANUP } from '../modules/local/cleanup/main'
include { UTILS_NFSCHEMA_PLUGIN } from '../subworkflows/nf-core/utils_nfschema_plugin/main'
include { BINNING } from '../subworkflows/local/binning'
include { BINNING_REFINEMENT } from '../subworkflows/local/binning_refinement'
include { BIN_QC } from '../subworkflows/local/bin_qc'
include { DOMAIN_CLASSIFICATION } from '../subworkflows/local/domain_classification'
include { BIN_ANNOTATION } from '../subworkflows/local/bin_annotation'
include { BIN_CLASSIFICATION } from '../subworkflows/local/bin_classification'

// Function to create input channel from CSV
def createCsvInputChannel(input_path) {
    return Channel
        .fromPath(input_path)
        .splitCsv(header:true, sep:',', strip:true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id ?: "sample_${System.currentTimeMillis()}"
            meta.run_id = row.run_id ?: 'default_run'
            meta.group = row.group ?: 'default_group'
            meta.single_end = false

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
                exit 1, "Missing or invalid paired-end read files for sample: ${meta.id}"
            }
            return tuple(meta, reads)
        }
}

workflow METAFLOW {
    main:
        // Parameter validation using UTILS_NFSCHEMA_PLUGIN
        UTILS_NFSCHEMA_PLUGIN (
            workflow,
            true,
            "${projectDir}/nextflow_schema.json"
        )

        // Preflight: normalize optional params + warnings
        def checkm2_db = params.checkm2_db?.toString()?.trim()
        if (!checkm2_db) {
            params.checkm2_db = null
            // Only show warning if user is actually using CheckM2 or Binette
            def needs_checkm2_db = (!params.enable_readbase && params.binqc_tool == 'checkm2' && !params.skip_binqc) || 
                                   (!params.enable_readbase && params.refine_tool == 'binette' && !params.skip_binning_refinement)
            if (needs_checkm2_db && !params.skip_binning && !params.skip_binning_bamabund) {
                log.warn "⚠️  Parameter --checkm2_db not provided; CheckM2 database will be downloaded automatically for CheckM2/Binette steps."
            }
        }

        // Create input channel
        input_ch = params.input_format == 'csv' ?
            createCsvInputChannel(params.input) :
            Channel.fromFilePairs("${params.input}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}")
                .map { sample_id, reads ->
                    def meta = [id:sample_id, single_end:false, run_id:'default_run', group:'default_group']
                    return tuple(meta, reads)
                }

        // Conditional preprocessing
        def cleaned_reads_source
        if (!params.skip_preprocess) {
            PREPROCESS(input_ch)
            cleaned_reads_source = PREPROCESS.out.cleaned_reads
        } else {
            // When skip_preprocess is true, use raw input reads directly
            cleaned_reads_source = input_ch
        }

        // Initialize empty channels with default values
        def read_based_versions_ch = Channel.empty()
        def read_based_results_ch = Channel.empty()
        def read_based_rgi_ch = Channel.empty()
        def assembly_versions_ch = Channel.empty()
        def assembly_megahit_contigs_ch = Channel.empty()
        def assembly_metaspades_contigs_ch = Channel.empty()
        def binning_bam_ch = Channel.empty()
        def binning_versions_ch = Channel.empty()
        def binning_bins_ch = Channel.empty()
        def binqc_versions_ch = Channel.empty()
        def binqc_summary_ch = Channel.empty()
        def binqc_quast_summary_ch = Channel.empty()
        def binqc_multiqc_ch = Channel.empty()
        def binannotation_versions_ch = Channel.empty()
        def binannotation_multiqc_ch = Channel.empty()
        def binclassification_versions_ch = Channel.empty()
        def binclassification_summary_ch = Channel.empty()

        // Use cleaned_reads_source for downstream workflows
        if (params.enable_readbase) {
            // Check if single-sample processing mode is enabled
            if (params.enable_singlesketch) {
                // Use READ_BASED_SINGLERUN for per-sample processing
                READ_BASED_SINGLERUN(cleaned_reads_source)
                read_based_versions_ch = READ_BASED_SINGLERUN.out.versions ?: Channel.empty()
                read_based_results_ch = READ_BASED_SINGLERUN.out.results ?: Channel.empty()
                read_based_rgi_ch = READ_BASED_SINGLERUN.out.rgi_results ?: Channel.empty()
            } else {
                // Use READ_BASED for batch processing (default)
                READ_BASED(cleaned_reads_source)
                read_based_versions_ch = READ_BASED.out.versions ?: Channel.empty()
                read_based_results_ch = READ_BASED.out.results ?: Channel.empty()
                read_based_rgi_ch = READ_BASED.out.rgi_results ?: Channel.empty()
            }
        } else {
            ASSEMBLY_BASED(cleaned_reads_source)
            assembly_versions_ch = ASSEMBLY_BASED.out.versions ?: Channel.empty()
            assembly_megahit_contigs_ch = ASSEMBLY_BASED.out.megahit_contigs ?: Channel.empty()
            assembly_metaspades_contigs_ch = ASSEMBLY_BASED.out.metaspades_contigs ?: Channel.empty()

            // Run BINNING_BAMABUND with assembly outputs
            def all_assemblies = assembly_megahit_contigs_ch.mix(assembly_metaspades_contigs_ch)
            
            if (!params.skip_binning_bamabund) {
                BINNING_BAMABUND(all_assemblies, cleaned_reads_source)
                binning_bam_ch = BINNING_BAMABUND.out.bam_bai ?: Channel.empty()
                binning_versions_ch = binning_versions_ch.mix(BINNING_BAMABUND.out.versions ?: Channel.empty())

                // Run BINNING with outputs from BINNING_BAMABUND
                if (!params.skip_binning && !params.skip_binning_bamabund) {
                    BINNING(all_assemblies, BINNING_BAMABUND.out.bam_bai)
                    binning_versions_ch = binning_versions_ch.mix(BINNING.out.versions ?: Channel.empty())

                    // Initialize channels for raw bins with metadata
                    def ch_raw_bins = BINNING.out.all_bins.map { meta, bin ->
                        def meta_new = meta.clone()
                        meta_new.refinement = 'unrefined'
                        [meta_new, bin]
                    }

                    // ============================================
                    // FIX: Run DOMAIN_CLASSIFICATION BEFORE refinement
                    // ============================================
                    if (params.bin_domain_classification) {
                        // Group raw bins by assembly for domain classification
                        ch_bins_for_classification = ch_raw_bins
                            .map { meta, bin ->
                                def meta_clean = meta.clone()
                                meta_clean.remove('bin_id')
                                def group_key = "${meta_clean.id}_${meta_clean.assembler}"
                                [group_key, meta_clean, bin]
                            }
                            .groupTuple(by: 0)
                            .map { key, metas, bins ->
                                [metas[0], bins.flatten()]
                            }

                        DOMAIN_CLASSIFICATION(all_assemblies, ch_bins_for_classification)
                        binning_versions_ch = binning_versions_ch.mix(DOMAIN_CLASSIFICATION.out.versions)

                        // Add domain metadata to raw bins
                        ch_domain_map = DOMAIN_CLASSIFICATION.out.classified_bins
                            .map { meta, bin ->
                                def bin_basename = bin.baseName
                                [bin_basename, meta.domain]
                            }

                        ch_raw_bins = ch_raw_bins
                            .map { meta, bin ->
                                def bin_basename = bin.baseName
                                [bin_basename, meta, bin]
                            }
                            .join(ch_domain_map, by: 0, remainder: true)
                            .map { bin_basename, meta, bin, domain ->
                                def meta_new = meta.clone()
                                meta_new.domain = domain ?: 'unclassified'
                                [meta_new, bin]
                            }
                    }

                    // ============================================
                    // Run BINNING_REFINEMENT (with domain already set)
                    // ============================================
                    def ch_refined_bins = Channel.empty()

                    if (!params.skip_binning_refinement) {

                        // Group bins for refinement (remove domain to group all together)
                        ch_bins_for_refinement = ch_raw_bins
                            .map { meta, bin ->
                                def meta_clean = meta.clone()
                                meta_clean.remove('domain')
                                meta_clean.remove('bin_id')
                                [meta_clean, bin]
                            }
                            .groupTuple()

                        BINNING_REFINEMENT(all_assemblies, ch_bins_for_refinement)
                        
                        // IMPORTANT: BINNING_REFINEMENT outputs are GROUPED bins, not individual bins
                        // We need to transpose them to match the raw bins channel format
                        ch_refined_bins = BINNING_REFINEMENT.out.refined_bins
                            .transpose()  // Add transpose here to get individual bins
                            .map { meta, bin ->
                                def meta_new = meta.clone()
                                meta_new.refinement = 'refined'
                                meta_new.bin_id = bin.baseName
                                // binrefine is already set by BINNING_REFINEMENT (DASTool or Binette)
                                [meta_new, bin]
                            }

                        // Re-add domain metadata to refined bins
                        if (params.bin_domain_classification) {
                            ch_refined_bins = ch_refined_bins
                                .map { meta, bin ->
                                    def bin_basename = bin.baseName
                                    [bin_basename, meta, bin]
                                }
                                .join(ch_domain_map, by: 0, remainder: true)
                                .map { bin_basename, meta, bin, domain ->
                                    def meta_new = meta.clone()
                                    meta_new.domain = domain ?: 'unclassified'
                                    [meta_new, bin]
                                }
                        }
                        
                        binning_versions_ch = binning_versions_ch.mix(BINNING_REFINEMENT.out.versions)
                    }

                    // Determine input for post-binning steps
                    if (!params.skip_binning_refinement) {
                        if (params.postbinning_input == 'raw_bins_only') {
                            binning_bins_ch = ch_raw_bins
                        } else if (params.postbinning_input == 'refined_bins_only') {
                            binning_bins_ch = ch_refined_bins
                        } else if (params.postbinning_input == 'both') {
                            // Both channels now have the same structure: [meta, single_bin]
                            binning_bins_ch = ch_raw_bins.mix(ch_refined_bins)
                        }
                    } else {
                        binning_bins_ch = ch_raw_bins
                    }

                    // Use multiMap to create separate channels for BIN_QC, BIN_ANNOTATION, and BIN_CLASSIFICATION
                    // This prevents channel consumption issues when processes run in parallel
                    ch_bins_split = binning_bins_ch
                        .multiMap { meta, bin ->
                            for_qc: [meta, bin]
                            for_annotation: [meta, bin]
                            for_classification: [meta, bin]
                        }

                    // Prepare input for BIN_QC - group bins
                    if (!params.skip_binqc && !params.skip_binning) {
                        ch_bins_for_qc = ch_bins_split.for_qc
                            .map { meta, bin ->
                                def meta_clean = meta.clone()
                                meta_clean.remove('bin_id')
                                meta_clean.remove('domain')
                                [meta_clean, bin]
                            }
                            .groupTuple()

                        BIN_QC(ch_bins_for_qc)
                        binqc_versions_ch = BIN_QC.out.versions ?: Channel.empty()
                        binqc_summary_ch = BIN_QC.out.qc_summary ?: Channel.empty()
                        binqc_quast_summary_ch = BIN_QC.out.quast_summary ?: Channel.empty()
                        binqc_multiqc_ch = BIN_QC.out.multiqc_files ?: Channel.empty()
                        binning_versions_ch = binning_versions_ch.mix(binqc_versions_ch)
                    }

                    // Run BIN_ANNOTATION in parallel with BIN_QC
                    if (!params.skip_bin_annotation && !params.skip_binning) {
                        ch_bins_for_annotation = ch_bins_split.for_annotation

                        BIN_ANNOTATION(ch_bins_for_annotation)
                        binannotation_versions_ch = BIN_ANNOTATION.out.versions ?: Channel.empty()
                        binannotation_multiqc_ch = BIN_ANNOTATION.out.multiqc_files ?: Channel.empty()
                        binning_versions_ch = binning_versions_ch.mix(binannotation_versions_ch)
                    }

                    // Run BIN_CLASSIFICATION in parallel with BIN_QC and BIN_ANNOTATION
                    if (!params.skip_bin_classification && !params.skip_binning) {
                        ch_bins_for_classification = ch_bins_split.for_classification

                        BIN_CLASSIFICATION(ch_bins_for_classification)
                        binclassification_versions_ch = BIN_CLASSIFICATION.out.versions ?: Channel.empty()
                        binclassification_summary_ch = BIN_CLASSIFICATION.out.genome_classification ?: Channel.empty()
                        binning_versions_ch = binning_versions_ch.mix(binclassification_versions_ch)
                    }
                }
            }
        }

        def combined_versions = Channel.empty()
        combined_versions = combined_versions.mix(assembly_versions_ch ?: Channel.empty())
        combined_versions = combined_versions.mix(read_based_versions_ch ?: Channel.empty())
        combined_versions = combined_versions.mix(binning_versions_ch ?: Channel.empty())

        // Optional cleanup
        if (params.cleanup) {
            def cleanup_sources = []
            if (params.enable_readbase) {
                cleanup_sources << read_based_results_ch
            } else {
                cleanup_sources << assembly_megahit_contigs_ch
                cleanup_sources << assembly_metaspades_contigs_ch
                cleanup_sources << binning_bam_ch
            }
            cleanup_sources = cleanup_sources.findAll { it }
            if (cleanup_sources) {
                def cleanup_trigger = Channel.merge(*cleanup_sources).collect()
                CLEANUP(cleanup_trigger)
            }
        }

    emit:
        versions = combined_versions.ifEmpty(null)
        results = read_based_results_ch.ifEmpty(null)
        read_based_rgi = read_based_rgi_ch.ifEmpty(null)
        assembly_megahit_contigs = assembly_megahit_contigs_ch.ifEmpty(null)
        assembly_metaspades_contigs = assembly_metaspades_contigs_ch.ifEmpty(null)
        binning_bam = binning_bam_ch.ifEmpty(null)
        binning_bins = binning_bins_ch.ifEmpty(null)
        binqc_summary = binqc_summary_ch.ifEmpty(null)
        binqc_quast_summary = binqc_quast_summary_ch.ifEmpty(null)
        binqc_multiqc = binqc_multiqc_ch.ifEmpty(null)
        binannotation_multiqc = binannotation_multiqc_ch.ifEmpty(null)
        binclassification_summary = binclassification_summary_ch.ifEmpty(null)
}
