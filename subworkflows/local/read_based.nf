#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules
include { SOURMASH_TAXANNOTATE as SOURMASH_TAXANNOTATE_META } from '../../modules/nf-core/sourmash/taxannotate/main'

// Import local modules
include { MERGE_PAIREDENDSEQS } from '../../modules/local/merge/main'
include { CREATE_BINS_CSV as CREATE_READS_CSV } from '../../modules/local/sourmash/createbinscsv/main'
include { SOURMASH_MANYSKETCH as SOURMASH_MANYSKETCH_META } from '../../modules/local/sourmash/manysketch/main'
include { EXTRACT_SOURMASH_SINGLESKETCHES as EXTRACT_SOURMASH_SINGLESKETCHES_META } from '../../modules/local/finalize/readbased/extract_sourmash_singlesketches'
include { SOURMASH_FASTMULTIGATHER as SOURMASH_FASTMULTIGATHER_META } from '../../modules/local/sourmash/fastmultigather/main'
include { SOURMASH_TAXMETAGENOME } from '../../modules/local/sourmash/taxmetagenome/main'
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
    // SOURMASH/YACHT Branch - New logic
    // ============================================================
    
    // Step 1: Merge paired-end reads
    MERGE_PAIREDENDSEQS(cleaned_reads_ch)
    versions_ch = versions_ch.mix(MERGE_PAIREDENDSEQS.out.versions)

    // Step 2: Create CSV manifest for merged sequences
    // Group all merged sequences together
    ch_merged_grouped = MERGE_PAIREDENDSEQS.out.merged_seqs
        .map { meta, merged_file ->
            def meta_batch = [id: 'batch']
            [meta_batch, merged_file]
        }
        .groupTuple()

    CREATE_READS_CSV(ch_merged_grouped)
    versions_ch = versions_ch.mix(CREATE_READS_CSV.out.versions.first())

    // Step 3: Run sourmash manysketch (new version)
    SOURMASH_MANYSKETCH_META(CREATE_READS_CSV.out.bins_csv)
    versions_ch = versions_ch.mix(SOURMASH_MANYSKETCH_META.out.versions.first())

    // Use multiMap to split sketch outputs for parallel processing
    ch_sketch_split = SOURMASH_MANYSKETCH_META.out.sketch_zip_file
        .multiMap { meta, sketch_zip ->
            for_fastmultigather: [meta, sketch_zip]
            for_yacht: [meta, sketch_zip]
            for_extract: [meta, sketch_zip]
        }

    // Step 4: Run sourmash fastmultigather (new version)
    SOURMASH_FASTMULTIGATHER_META(
        ch_sketch_split.for_fastmultigather,
        file(params.sourmash_database, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_FASTMULTIGATHER_META.out.versions.first())

    // Step 5: Run YACHT (parallel with step 6-8)
    // Extract zip files directory for YACHT
    ch_yacht_input = ch_sketch_split.for_yacht
        .map { meta, sketch_zip ->
            // YACHT needs directory containing individual zip files
            // We'll use EXTRACT_SOURMASH_SINGLESKETCHES output
            [meta, sketch_zip]
        }

    // Step 9: Extract single sketches (parallel branch - for publishing only)
    EXTRACT_SOURMASH_SINGLESKETCHES_META(
        ch_sketch_split.for_extract,
        params.sourmash_ksize
    )
    versions_ch = versions_ch.mix(EXTRACT_SOURMASH_SINGLESKETCHES_META.out.versions.first())

    // Use the extracted zip files directory for YACHT
    YACHT_RUN(
        EXTRACT_SOURMASH_SINGLESKETCHES_META.out.manysketch_dir.map { meta, dir -> 
            [meta, file("${dir}/zip_files")]
        },
        file(params.yacht_database, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(YACHT_RUN.out.versions.first())

    // Step 6: Run sourmash tax annotate
    SOURMASH_TAXANNOTATE_META(
        SOURMASH_FASTMULTIGATHER_META.out.gather_csv,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_TAXANNOTATE_META.out.versions.first())

    // Decompress the .with-lineages.csv.gz file
    ch_taxannotate_decompressed = SOURMASH_TAXANNOTATE_META.out.result
        .map { meta, csv_gz ->
            // The process will output .csv.gz, we need to decompress it
            [meta, csv_gz]
        }

    // Step 7: Process results (uses decompressed taxannotate and YACHT results)
    PROCESS_READBASED_RESULTS(
        ch_taxannotate_decompressed,
        YACHT_RUN.out.yacht_xlsx
    )
    versions_ch = versions_ch.mix(PROCESS_READBASED_RESULTS.out.versions.first())

    // Step 8: Run sourmash tax metagenome
    SOURMASH_TAXMETAGENOME(
        SOURMASH_TAXANNOTATE_META.out.result,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    versions_ch = versions_ch.mix(SOURMASH_TAXMETAGENOME.out.versions.first())

    emit:
    versions = versions_ch.ifEmpty(null)
    results = PROCESS_READBASED_RESULTS.out.final_results
    rgi_results = ch_rgi_results
    gather_csv = SOURMASH_FASTMULTIGATHER_META.out.gather_csv
    taxannotate = SOURMASH_TAXANNOTATE_META.out.result
    metagenome_classification = SOURMASH_TAXMETAGENOME.out.genome_classification
    single_sketches = EXTRACT_SOURMASH_SINGLESKETCHES_META.out.sig_zip_files
}