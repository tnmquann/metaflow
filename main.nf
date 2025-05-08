#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { MERGE_PAIREDENDSEQS } from './modules/common/local/merge/main'
include { SOURMASH_MANYSKETCH } from './modules/read_based/local/sourmash/manysketch'
include { SOURMASH_FASTMULTIGATHER } from './modules/read_based/local/sourmash/fastmultigather'
include { YACHT_RUN } from './modules/read_based/local/yacht/run'
include { PROCESS_READBASED_RESULTS } from './modules/common/local/results/process_readbased'
include { CLEANUP } from './modules/common/local/cleanup/main'

// Function to create input channel based on CSV
def create_input_channel_from_csv(csv_file) {
    Channel
        .fromPath(csv_file)
        .splitCsv(header:true)
        .map { row -> 
            def meta = [
                id: row.sample_id,
                run_id: row.run_id,
                group: row.group
            ]
            def reads = []
            if (row.short_reads_1 && row.short_reads_2) {
                reads = [
                    file(row.short_reads_1),
                    file(row.short_reads_2)
                ]
            }
            if (row.long_reads) {
                reads << file(row.long_reads)
            }
            return tuple(meta, reads)
        }
}

// Function to create input channel from directory
def create_input_channel_from_dir(dir_path) {
    Channel
        .fromFilePairs("${dir_path}/*_{1,2}.{fastq,fq}.gz")
        .map { id, files -> 
            def meta = [
                id: id,
                run_id: '0',
                group: '0'
            ]
            return tuple(meta, files)
        }
}

workflow {
    main:
        // Create input channel based on input format
        ch_input = params.input_format == 'csv' && params.input ? 
            create_input_channel_from_csv(params.input) :
            create_input_channel_from_dir(params.trimmed_fastq)

        // Run MERGE_PAIREDENDSEQS on each sample
        MERGE_PAIREDENDSEQS(ch_input)

        // Group merged files by batch for SOURMASH_MANYSKETCH
        ch_merged_batch = MERGE_PAIREDENDSEQS.out.merged_seqs
            .map { meta, file -> file }
            .collect()
            .map { files ->
                def content = ["name,genome_filename,protein_filename"]
                files.each { file ->
                    def name = file.name.replace('.fastq.gz', '')
                    content << "${name},${file},"
                }
                def csv = file("${workflow.workDir}/manysketch_input.csv")
                csv.text = content.join('\n')
                return csv
            }

        // Continue with existing workflow
        SOURMASH_MANYSKETCH(ch_merged_batch)

        ch_sourmash_sketches = SOURMASH_MANYSKETCH.out.sketch_zip_file
            .map { sketch_zip -> 
                [ [id: 'batch'], sketch_zip ] 
            }

        ch_sourmash_signatures = SOURMASH_MANYSKETCH.out.zip_files_dir
            .map { zip_dir -> 
                [ [id: 'batch'], zip_dir ] 
            }

        SOURMASH_FASTMULTIGATHER(
            ch_sourmash_sketches,
            file(params.sourmash_database)
        )

        YACHT_RUN(
            ch_sourmash_signatures,
            file(params.yacht_database)
        )

        PROCESS_READBASED_RESULTS(
            SOURMASH_FASTMULTIGATHER.out.gather_csv,
            YACHT_RUN.out.yacht_xlsx
        )

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