#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESS } from '../subworkflows/local/preprocess'
include { READ_BASED } from '../subworkflows/local/read_based'
include { CLEANUP } from '../modules/local/cleanup/main'
include { UTILS_NFSCHEMA_PLUGIN } from '../subworkflows/nf-core/utils_nfschema_plugin/main'

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
        READ_BASED(PREPROCESS.out.cleaned_reads)

        // Optional cleanup
        if (params.cleanup) {
            CLEANUP(READ_BASED.out.results.collect())
        }

    emit:
        versions = READ_BASED.out.versions // Emit versions for tracking
        results = READ_BASED.out.results  // Emit final results
}