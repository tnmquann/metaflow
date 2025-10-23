#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESS } from '../subworkflows/local/preprocess'
include { READ_BASED } from '../subworkflows/local/read_based'
include { ASSEMBLY_BASED } from '../subworkflows/local/assembly_based'
include { BINNING_BAMABUND } from '../subworkflows/local/binning_bamabund'
include { CLEANUP } from '../modules/local/cleanup/main'
include { UTILS_NFSCHEMA_PLUGIN } from '../subworkflows/nf-core/utils_nfschema_plugin/main'
include { BINNING } from '../subworkflows/local/binning'
include { BINNING_REFINEMENT } from '../subworkflows/local/binning_refinement'

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
            workflow,           // Pass workflow object
            true,              // Validate parameters
            "${projectDir}/nextflow_schema.json" // Schema path
        )

        // Create input channel
        input_ch = params.input_format == 'csv' ?
            createCsvInputChannel(params.input) :
            Channel.fromFilePairs("${params.input}/*_{1,2}*.fastq.gz")
                .map { sample_id, reads ->
                    def meta = [id:sample_id, single_end:false, run_id:'default_run', group:'default_group']
                    return tuple(meta, reads)
                }

        // Run subworkflows
        PREPROCESS(input_ch)

        def cleaned_reads_source = PREPROCESS.out.cleaned_reads
        def read_based_versions_ch = Channel.empty()
        def read_based_results_ch = Channel.empty()
        def read_based_rgi_ch = Channel.empty()

        def assembly_versions_ch = Channel.empty()
        def assembly_megahit_contigs_ch = Channel.empty()
        def assembly_metaspades_contigs_ch = Channel.empty()
        def binning_bam_ch = Channel.empty()
        def binning_versions_ch = Channel.empty()
        def binning_bins_ch = Channel.empty()

        if (params.enable_readbase) {
            READ_BASED(cleaned_reads_source)
            read_based_versions_ch = READ_BASED.out.versions ?: Channel.empty()
            read_based_results_ch = READ_BASED.out.results ?: Channel.empty()
            read_based_rgi_ch = READ_BASED.out.rgi_results ?: Channel.empty()
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
                binning_versions_ch = BINNING_BAMABUND.out.versions ?: Channel.empty()

                // Run BINNING with outputs from BINNING_BAMABUND
                if (!params.skip_binning && !params.skip_binning_bamabund) {
                    BINNING(all_assemblies, BINNING_BAMABUND.out.bam_bai)
                    binning_bins_ch = BINNING.out.all_bins ?: Channel.empty()
                    binning_versions_ch = binning_versions_ch.mix(BINNING.out.versions ?: Channel.empty())

                    // Add binning refinement
                    if (!params.skip_binning_refinement && params.refine_tool == 'dastool') {
                        BINNING_REFINEMENT(all_assemblies, BINNING.out.all_bins)
                        
                        // Choose which bins to use for downstream analysis
                        if (params.postbinning_input == 'refined_bins_only') {
                            binning_bins_ch = BINNING_REFINEMENT.out.refined_bins
                        } else if (params.postbinning_input == 'both') {
                            binning_bins_ch = BINNING.out.all_bins.mix(BINNING_REFINEMENT.out.refined_bins)
                        }
                        // If 'raw_bins_only', keep original binning_bins_ch
                        
                        binning_versions_ch = binning_versions_ch.mix(BINNING_REFINEMENT.out.versions)
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
}
