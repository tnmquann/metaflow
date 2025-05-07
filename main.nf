#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { MERGE_SEQUENCES } from './modules/merge_sequences'
include { SOURMASH_MANYSKETCH } from './modules/sourmash/manysketch'
include { SOURMASH_FASTMULTIGATHER } from './modules/sourmash/fastmultigather'
include { YACHT } from './modules/yacht/yacht'
include { PROCESS_RESULTS } from './modules/process_results'
include { CLEANUP } from './modules/cleanup'

// Main workflow
workflow {
    main:
        // Create channel for input directory
        ch_input_dir = Channel.fromPath(params.trimmed_fastq)

        // Run MERGE_SEQUENCES on the input directory
        MERGE_SEQUENCES(ch_input_dir)

        // Create the manysketch.csv file content
        ch_manysketch_csv = MERGE_SEQUENCES.out.merged_seqs
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
        SOURMASH_MANYSKETCH(ch_manysketch_csv)

        // Create metadata for batch processing - using the zip file
        ch_sketch_meta = SOURMASH_MANYSKETCH.out.sketch_zip_file
            .map { sketch_zip -> 
                [ [id: 'batch'], sketch_zip ] 
            }

        // Using the zip_files_dir output
        ch_zip_files_meta = SOURMASH_MANYSKETCH.out.zip_files_dir
            .map { zip_dir -> 
                [ [id: 'batch'], zip_dir ] 
            }

        // Run SOURMASH_FASTMULTIGATHER with metadata
        SOURMASH_FASTMULTIGATHER(
            ch_sketch_meta,
            file(params.sourmash_database)
        )

        // Run YACHT with metadata
        YACHT(
            ch_zip_files_meta,
            file(params.yacht_database)
        )

        // Process results
        PROCESS_RESULTS(
            SOURMASH_FASTMULTIGATHER.out.gather_csv,
            YACHT.out.yacht_xlsx
        )

        // Run cleanup if enabled
        if (params.cleanup) {
            CLEANUP(PROCESS_RESULTS.out.final_results.collect())
        }
}

// Update subworkflows accordingly
workflow SOURMASH_ONLY {
    take:
        input_dir    // Directory path

    main:
        MERGE_SEQUENCES(input_dir)

        ch_manysketch_csv = MERGE_SEQUENCES.out.merged_seqs
            .map { merged_seq_dir ->
                def content = []
                content << "name,genome_filename,protein_filename"
                file(merged_seq_dir).listFiles().each { file ->
                    if (file.name.endsWith('.fastq.gz')) {
                        def name = file.name.replace('.fastq.gz', '')
                        content << "${name},${file.toString()},"
                    }
                }
                return content.join('\n')
            }
            .map { content ->
                def csv = file("${workflow.workDir}/manysketch_input_subworkflow.csv")
                csv.text = content
                return csv
            }

        SOURMASH_MANYSKETCH(ch_manysketch_csv)

        // Create metadata for batch processing
        ch_sketch_meta = SOURMASH_MANYSKETCH.out.sketch_zip_file
            .map { sketch_zip -> 
                [ [id: 'batch'], sketch_zip ] 
            }

        SOURMASH_FASTMULTIGATHER(
            ch_sketch_meta,
            file(params.sourmash_database)
        )

    emit:
        results = SOURMASH_FASTMULTIGATHER.out.gather_csv
        manysketch_sketch_zip = SOURMASH_MANYSKETCH.out.sketch_zip_file
        manysketch_zip_files_dir = SOURMASH_MANYSKETCH.out.zip_files_dir
}

workflow YACHT_ONLY {
    take:
        input_dir    // Directory path

    main:
        sourmash_outputs = SOURMASH_ONLY(input_dir)

        ch_zip_files_meta = SOURMASH_MANYSKETCH.out.zip_files_dir
            .map { zip_dir -> 
                [ [id: 'batch'], zip_dir ] 
            }

        YACHT(
            ch_zip_files_meta,
            file(params.yacht_database)
        )

    emit:
        results = YACHT.out.yacht_xlsx
}