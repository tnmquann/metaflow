#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { MERGE_PAIREDENDSEQS } from './modules/common/local/merge/main'
include { SOURMASH_MANYSKETCH } from './modules/read_based/local/sourmash/manysketch'
include { SOURMASH_FASTMULTIGATHER } from './modules/read_based/local/sourmash/fastmultigather'
include { YACHT_RUN } from './modules/read_based/local/yacht/run'
include { PROCESS_READBASED_RESULTS } from './modules/common/local/results/process_readbased'
include { CLEANUP } from './modules/common/local/cleanup/main'

// Main workflow
workflow {
    main:
        // Create channel for input directory
        ch_input_fastq = Channel.fromPath(params.trimmed_fastq)

        // Run MERGE_PAIREDENDSEQS on the input directory
        MERGE_PAIREDENDSEQS(ch_input_fastq)

        // Create the manysketch.csv file content
        ch_sourmash_manifest = MERGE_PAIREDENDSEQS.out.merged_seqs
            .map { merged_seq_dir ->
                def content = []
                content << "name,genome_filename,protein_filename"
                
                // Process each merged file
                file(merged_seq_dir).listFiles().each { file ->
                    if (file.name.endsWith('.fastq.gz')) {
                        def name = file.name.replace('.fastq.gz', '')
                        content << "${name},${file.toString()},"
                    }
                }
                return content.join('\n')
            }
            .map { content ->
                def csv = file("${workflow.workDir}/manysketch_input.csv")
                csv.text = content
                return csv
            }

        // Run SOURMASH_MANYSKETCH with the CSV file
        SOURMASH_MANYSKETCH(ch_sourmash_manifest)

        // Create metadata for batch processing - using the zip file
        ch_sourmash_sketches = SOURMASH_MANYSKETCH.out.sketch_zip_file
            .map { sketch_zip -> 
                [ [id: 'batch'], sketch_zip ] 
            }

        // Using the zip_files_dir output
        ch_sourmash_signatures = SOURMASH_MANYSKETCH.out.zip_files_dir
            .map { zip_dir -> 
                [ [id: 'batch'], zip_dir ] 
            }

        // Run SOURMASH_FASTMULTIGATHER with metadata
        SOURMASH_FASTMULTIGATHER(
            ch_sourmash_sketches,
            file(params.sourmash_database)
        )

        // Run YACHT_RUN with metadata
        YACHT_RUN(
            ch_sourmash_signatures,
            file(params.yacht_database)
        )

        // Process results
        PROCESS_READBASED_RESULTS(
            SOURMASH_FASTMULTIGATHER.out.gather_csv,
            YACHT_RUN.out.yacht_xlsx
        )

        // Run cleanup if enabled
        if (params.cleanup) {
            CLEANUP(PROCESS_READBASED_RESULTS.out.final_results.collect())
        }
}

// Update subworkflow names
workflow PROFILING_SOURMASHONLY {
    take:
        ch_input_fastq    // More descriptive parameter name

    main:
        MERGE_PAIREDENDSEQS(ch_input_fastq)

        ch_sourmash_manifest = MERGE_PAIREDENDSEQS.out.merged_seqs
            .map { createSourmashManifest(it) }

        SOURMASH_MANYSKETCH(ch_sourmash_manifest)

        ch_sourmash_sketches = SOURMASH_MANYSKETCH.out.sketch_zip_file
            .map { sketch_zip -> 
                [ [id: 'batch'], sketch_zip ] 
            }

        SOURMASH_FASTMULTIGATHER(
            ch_sourmash_sketches,
            file(params.sourmash_database)
        )

    emit:
        results = SOURMASH_FASTMULTIGATHER.out.gather_csv
        sketches = SOURMASH_MANYSKETCH.out.sketch_zip_file
        signatures = SOURMASH_MANYSKETCH.out.zip_files_dir
}

workflow PROFILING_YACHTONLY {
    take:
        ch_input_fastq    // More descriptive parameter name

    main:
        sourmash_results = PROFILING_SOURMASHONLY(ch_input_fastq)

        ch_sourmash_signatures = sourmash_results.signatures
            .map { sig_dir -> 
                [ [id: 'batch'], sig_dir ] 
            }

        YACHT_RUN(
            ch_sourmash_signatures,
            file(params.yacht_database)
        )

    emit:
        results = YACHT_RUN.out.yacht_xlsx
}

// Helper function for creating Sourmash manifest
def createSourmashManifest(merged_seq_dir) {
    def content = []
    content << "name,genome_filename,protein_filename"
    file(merged_seq_dir).listFiles().each { file ->
        if (file.name.endsWith('.fastq.gz')) {
            def name = file.name.replace('.fastq.gz', '')
            content << "${name},${file.toString()},"
        }
    }
    def csv = file("${workflow.workDir}/manysketch_input_subworkflow.csv")
    csv.text = content.join('\n')
    return csv
}