#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIN_CLASSIFICATION SUBWORKFLOW
    Taxonomic classification of bins using sourmash
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import nf-core modules
include { SOURMASH_TAXANNOTATE }            from '../../modules/nf-core/sourmash/taxannotate/main'

// Import local modules
include { CREATE_BINS_CSV }                         from '../../modules/local/sourmash/createbinscsv/main'
include { SOURMASH_MANYSKETCH as SOURMASH_MANYSKETCH_BINS }     from '../../modules/local/sourmash/manysketch2/main'
include { EXTRACT_SOURMASH_SINGLESKETCHES }         from '../../modules/local/finalize/assemblybased/extract_sourmash_singlesketches'
include { SOURMASH_FASTMULTIGATHER as SOURMASH_FASTMULTIGATHER_BINS } from '../../modules/local/sourmash/fastmultigather2/main'
include { COMBINE_FASTMULTIGATHER_RESULTS } from '../../modules/local/finalize/assemblybased/combine_fastmultigather_results'
include { SOURMASH_TAXGENOME }              from '../../modules/local/sourmash/taxgenome/main'

workflow BIN_CLASSIFICATION {
    take:
    bins               // channel: [ val(meta), path(bin.fa) ] - individual bins from BINNING/BINNING_REFINEMENT

    main:
    ch_versions = Channel.empty()

    // Validate required parameters
    if (!params.sourmash_database) {
        error "ERROR: Please provide --sourmash_database parameter for bin classification"
    }
    if (!params.sourmash_taxonomy_csv) {
        error "ERROR: Please provide --sourmash_taxonomy_csv parameter for bin classification"
    }

    // ============================================
    // Group bins by assembly/binner/refinement for batch processing
    // ============================================
    ch_bins_grouped = bins
        .map { meta, bin ->
            def meta_clean = meta.clone()
            meta_clean.remove('bin_id')
            meta_clean.remove('domain')
            [meta_clean, bin]
        }
        .groupTuple()
        .map { meta, bin_files ->
            [meta, bin_files.flatten()]
        }

    // Step 1: Create CSV file listing all bins
    CREATE_BINS_CSV(ch_bins_grouped)
    ch_versions = ch_versions.mix(CREATE_BINS_CSV.out.versions)

    // Step 2: Generate signatures using manysketch WITH metadata
    SOURMASH_MANYSKETCH_BINS(CREATE_BINS_CSV.out.bins_csv)
    ch_versions = ch_versions.mix(SOURMASH_MANYSKETCH_BINS.out.versions)

    // Step 2b: Extract single sketches from the batch sketch zip
    EXTRACT_SOURMASH_SINGLESKETCHES(
        SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file.map { meta, zip -> zip },
        params.sourmash_ksize_bins
    )
    ch_versions = ch_versions.mix(EXTRACT_SOURMASH_SINGLESKETCHES.out.versions)

    // Step 3: Prepare sketch with metadata for fastmultigather
    ch_sketch_for_gather = SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file

    // Step 4: Handle multiple databases (space-separated)
    if (params.sourmash_database.contains(' ')) {
        db_files = params.sourmash_database.split(' ')
        ch_sourmash_database = Channel.fromList(db_files.collect { it.trim() })
            .map { db_path -> file(db_path, checkIfExists: true) }
    } else {
        ch_sourmash_database = Channel.fromPath(params.sourmash_database, checkIfExists: true)
    }

    // Step 5: Run fastmultigather for each database
    ch_sketch_database_combinations = ch_sketch_for_gather
        .combine(ch_sourmash_database)
        .map { meta, sketch, database -> 
            [meta, sketch, database]
        }

    SOURMASH_FASTMULTIGATHER_BINS(
        ch_sketch_database_combinations.map { meta, sketch, db -> [meta, sketch] },
        ch_sketch_database_combinations.map { meta, sketch, db -> db }
    )
    ch_versions = ch_versions.mix(SOURMASH_FASTMULTIGATHER_BINS.out.versions)

    // ============================================
    // Step 6: COMBINE results BY meta.id
    // Group all CSV files from different assembler/binner/binrefine by meta.id
    // ============================================
    ch_combined_results = SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv
        .map { meta, csv_file ->
            // Group by meta.id only, removing assembler/binner/binrefine distinctions
            [meta.id, meta, csv_file]
        }
        .groupTuple(by: 0)
        .map { id, metas, csv_files ->
            // Use the first meta as template, update with combined info
            def meta_combined = metas[0].clone()
            meta_combined.id = id
            // Remove specific assembler/binner/binrefine to indicate this is combined
            meta_combined.remove('assembler')
            meta_combined.remove('binner')
            meta_combined.remove('binrefine')
            [meta_combined, csv_files.flatten()]
        }

    COMBINE_FASTMULTIGATHER_RESULTS(ch_combined_results)
    ch_versions = ch_versions.mix(COMBINE_FASTMULTIGATHER_RESULTS.out.versions)
    ch_final_gather_csv = COMBINE_FASTMULTIGATHER_RESULTS.out.gather_csv

    // Step 7: Taxonomy annotation
    ch_sourmash_taxonomy_csv = Channel.fromPath(params.sourmash_taxonomy_csv, checkIfExists: true)

    SOURMASH_TAXANNOTATE(
        ch_final_gather_csv,
        ch_sourmash_taxonomy_csv
    )
    ch_versions = ch_versions.mix(SOURMASH_TAXANNOTATE.out.versions)

    // Step 8: Generate genome classification
    SOURMASH_TAXGENOME(
        SOURMASH_TAXANNOTATE.out.result.map { meta, csv -> 
            [meta, csv]
        },
        ch_sourmash_taxonomy_csv
    )
    ch_versions = ch_versions.mix(SOURMASH_TAXGENOME.out.versions)

    emit:
    bins_csv            = CREATE_BINS_CSV.out.bins_csv
    sketch_zip          = SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file
    single_sketches_dir = EXTRACT_SOURMASH_SINGLESKETCHES.out.manysketch_dir
    single_sketches_zip = EXTRACT_SOURMASH_SINGLESKETCHES.out.zip_files_dir
    gather_csv          = ch_final_gather_csv
    annotated_csv       = SOURMASH_TAXANNOTATE.out.result
    genome_classification = SOURMASH_TAXGENOME.out.genome_classification
    multiqc_files       = Channel.empty()
    versions            = ch_versions
}
