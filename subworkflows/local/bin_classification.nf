#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules
include { SOURMASH_TAXANNOTATE }            from '../../modules/nf-core/sourmash/taxannotate/main'

// Import local modules
include { CREATE_BINS_CSV }                         from '../../modules/local/sourmash/createbinscsv/main'
include { SOURMASH_MANYSKETCH as SOURMASH_MANYSKETCH_BINS }     from '../../modules/local/sourmash/manysketch/main'
include { EXTRACT_SOURMASH_SINGLESKETCHES }         from '../../modules/local/finalize/readbased/extract_sourmash_singlesketches'
include { SOURMASH_FASTMULTIGATHER as SOURMASH_FASTMULTIGATHER_BINS } from '../../modules/local/sourmash/fastmultigather/main'
include { COMBINE_FASTMULTIGATHER_RESULTS } from '../../modules/local/finalize/readbased/combine_fastmultigather_results'
include { COMBINE_FASTMULTIGATHER_RESULTS as COMBINE_TAXANNOTATE_RESULTS} from '../../modules/local/finalize/readbased/combine_fastmultigather_results'
include { COMBINE_FASTMULTIGATHER_RESULTS as COMBINE_TAXGENOME_RESULTS} from '../../modules/local/finalize/readbased/combine_fastmultigather_results'
include { SOURMASH_TAXGENOME }              from '../../modules/local/sourmash/taxgenome/main'

workflow BIN_CLASSIFICATION {
    take:
    bins_ch  // channel: [ val(meta), path(bin) ] - individual bins from BINNING/BINNING_REFINEMENT

    main:
    ch_versions = Channel.empty()

    // Validate required parameters
    if (!params.sourmash_database) {
        error "sourmash_database parameter is required for bin classification"
    }
    if (!params.sourmash_taxonomy_csv) {
        error "sourmash_taxonomy_csv parameter is required for bin classification"
    }

    // Step 1: Group bins by assembly/binner/binrefine/sample to create CSV manifests
    ch_bins_grouped = bins_ch
        .map { meta, bin ->
            def meta_clean = meta.clone()
            meta_clean.remove('bin_id')
            meta_clean.remove('domain')
            [meta_clean, bin]
        }
        .groupTuple()

    // Create CSV files for each group
    CREATE_BINS_CSV(ch_bins_grouped)
    ch_versions = ch_versions.mix(CREATE_BINS_CSV.out.versions.first())

    // Step 2: Run sourmash manysketch on CSV files
    SOURMASH_MANYSKETCH_BINS(CREATE_BINS_CSV.out.bins_csv)
    ch_versions = ch_versions.mix(SOURMASH_MANYSKETCH_BINS.out.versions.first())

    // Use multiMap to split the sketch outputs for parallel processing
    ch_sketch_split = SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file
        .multiMap { meta, sketch_zip ->
            for_fastmultigather: [meta, sketch_zip]
            for_extract: [meta, sketch_zip]
        }

    // Step 3: Run fastmultigather
    // Output: ${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}_sourmash_gather.csv
    SOURMASH_FASTMULTIGATHER_BINS(
        ch_sketch_split.for_fastmultigather,
        file(params.sourmash_database, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SOURMASH_FASTMULTIGATHER_BINS.out.versions.first())

    // Step 9: Extract single sketches (parallel branch - only for publishing)
    EXTRACT_SOURMASH_SINGLESKETCHES(
        ch_sketch_split.for_extract,
        params.sourmash_ksize
    )
    ch_versions = ch_versions.mix(EXTRACT_SOURMASH_SINGLESKETCHES.out.versions.first())

    // Step 4: Run sourmash tax annotate
    // Input: gather_csv from step 3
    // Output: ${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}_sourmash_gather.with-lineages.csv.gz
    SOURMASH_TAXANNOTATE(
        SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SOURMASH_TAXANNOTATE.out.versions.first())

    // Step 5: Run sourmash tax genome
    // Input: TAXANNOTATE output (with-lineages.csv.gz)
    // Output: ${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}_genome_classification.csv
    SOURMASH_TAXGENOME(
        SOURMASH_TAXANNOTATE.out.result,
        file(params.sourmash_taxonomy_csv, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SOURMASH_TAXGENOME.out.versions.first())

    // Step 6: Combine fastmultigather results by sample ID
    ch_gather_for_combine = SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv
        .map { meta, csv ->
            def meta_combine = [:]
            meta_combine.id = meta.id
            meta_combine.file_type = 'gather'
            [meta_combine.id, meta_combine, csv]
        }
        .groupTuple(by: 0)
        .map { id, metas, csvs ->
            [metas[0], csvs.flatten()]
        }

    COMBINE_FASTMULTIGATHER_RESULTS(ch_gather_for_combine)
    ch_versions = ch_versions.mix(COMBINE_FASTMULTIGATHER_RESULTS.out.versions.first())

    // Step 7: Combine taxannotate results by sample ID
    // Note: Decompress .gz files before combining
    ch_taxannotate_for_combine = SOURMASH_TAXANNOTATE.out.result
        .map { meta, csv_gz ->
            def meta_combine = [:]
            meta_combine.id = meta.id
            meta_combine.file_type = 'taxannotate'
            [meta_combine.id, meta_combine, csv_gz]
        }
        .groupTuple(by: 0)
        .map { id, metas, csvs ->
            [metas[0], csvs.flatten()]
        }

    COMBINE_TAXANNOTATE_RESULTS(ch_taxannotate_for_combine)
    ch_versions = ch_versions.mix(COMBINE_TAXANNOTATE_RESULTS.out.versions.first())

    // Step 8: Combine taxgenome results by sample ID
    ch_taxgenome_for_combine = SOURMASH_TAXGENOME.out.genome_classification
        .map { meta, csv ->
            def meta_combine = [:]
            meta_combine.id = meta.id
            meta_combine.file_type = 'taxgenome'
            [meta_combine.id, meta_combine, csv]
        }
        .groupTuple(by: 0)
        .map { id, metas, csvs ->
            [metas[0], csvs.flatten()]
        }

    COMBINE_TAXGENOME_RESULTS(ch_taxgenome_for_combine)
    ch_versions = ch_versions.mix(COMBINE_TAXGENOME_RESULTS.out.versions.first())

    emit:
    bins_csv             = CREATE_BINS_CSV.out.bins_csv
    sketches             = SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file
    single_sketches      = EXTRACT_SOURMASH_SINGLESKETCHES.out.sig_zip_files
    gather_csv           = SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv
    taxannotate          = SOURMASH_TAXANNOTATE.out.result
    genome_classification = SOURMASH_TAXGENOME.out.genome_classification
    combined_gather      = COMBINE_FASTMULTIGATHER_RESULTS.out.combined_csv
    combined_taxannotate = COMBINE_TAXANNOTATE_RESULTS.out.combined_csv
    combined_taxgenome   = COMBINE_TAXGENOME_RESULTS.out.combined_csv
    versions             = ch_versions
}
