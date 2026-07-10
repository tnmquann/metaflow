#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESS } from '../subworkflows/local/preprocess/main'
include { READ_BASED } from '../subworkflows/local/read_based/main'
include { READ_BASED_SINGLERUN } from '../subworkflows/local/read_based_singlerun/main'
include { ASSEMBLY_BASED } from '../subworkflows/local/assembly_based/main'
include { BINNING_BAMABUND } from '../subworkflows/local/binning_bamabund/main'
include { CLEANUP } from '../modules/local/cleanup/main'
include { UTILS_NFSCHEMA_PLUGIN } from '../subworkflows/nf-core/utils_nfschema_plugin/main'
include { BINNING } from '../subworkflows/local/binning/main'
include { BINNING_REFINEMENT } from '../subworkflows/local/binning_refinement/main'
include { BIN_QC } from '../subworkflows/local/bin_qc/main'
include {
    DOMAIN_CLASSIFICATION as DOMAIN_CLASSIFICATION_RAW
    DOMAIN_CLASSIFICATION as DOMAIN_CLASSIFICATION_REFINED
} from '../subworkflows/local/domain_classification/main'
include { BIN_ANNOTATION } from '../subworkflows/local/bin_annotation/main'
include { BIN_CLASSIFICATION } from '../subworkflows/local/bin_classification/main'
include { GUNZIP as GUNZIP_ASSEMBLIES } from '../modules/nf-core/gunzip/main'

// Import modules for database setup
include { CHECKM2_DATABASEDOWNLOAD } from '../modules/nf-core/checkm2/databasedownload/main'

// Function to create input channel from CSV
def createCsvInputChannel(input_path) {
    return channel
        .fromPath(input_path)
        .splitCsv(header:true, sep:',', strip:true)
        .map { row ->
            def meta = [:]
            def sample_id = row.sample_id?.toString()?.trim()
            if (!sample_id) {
                exit 1, 'Missing required sample_id in CSV input.'
            }
            meta.id = sample_id
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

// Function to create assembly input channel from CSV
def createAssemblyInputChannel(input_path) {
    return channel
        .fromPath(input_path)
        .splitCsv(header:true, sep:',', strip:true)
        .map { row ->
            def meta = [:]
            def sample_id = row.sample_id?.toString()?.trim()
            if (!sample_id) {
                exit 1, 'Missing required sample_id in assembly CSV input.'
            }
            meta.id = sample_id
            meta.run_id = row.run_id ?: 'default_run'
            meta.group = row.group ?: 'default_group'
            
            def assembler = row.assembler?.trim()
            if (!assembler) {
                exit 1, "Missing assembler for sample ${meta.id}. Please specify the assembler used (e.g., MEGAHIT, SPAdes, IDBA-UD, etc.)"
            }
            meta.assembler = assembler
            
            def fasta_path = row.fasta?.trim()
            if (!fasta_path) {
                exit 1, "Missing fasta file path for sample ${meta.id}"
            }
            
            def fasta_file = file(fasta_path)
            if (!fasta_file.exists()) {
                exit 1, "Fasta file does not exist for sample ${meta.id}: ${fasta_path}"
            }
            
            return tuple(meta, fasta_file)
        }
}


workflow METAFLOW {
    main:
        schema_path = projectDir.name == 'workflow' ? '../nextflow_schema.json' : 'nextflow_schema.json'

        // Normalize nullable option strings before nf-schema validation.
        def readbased_postprocess = params.readbased_postprocess?.toString()?.trim()
        params.readbased_postprocess = readbased_postprocess ?: null
        params.postprocess_options = params.postprocess_options?.toString()?.trim() ?: ''

        // Parameter validation using UTILS_NFSCHEMA_PLUGIN
        UTILS_NFSCHEMA_PLUGIN (
            workflow,
            !workflow.preview,
            schema_path,
            params.help ?: false,
            params.help_full ?: false,
            params.show_hidden ?: false,
            '',
            '',
            '',
            true
        )

        // Preflight: validate optional read-based post-processing dependencies.
        if (params.readbased_postprocess) {
            def valid_postprocess_modes = ['all', 'phyloseq', 'taxburst', 'rgi_bwt']
            if (!params.enable_readbase) {
                exit 1, 'ERROR: --readbased_postprocess can only be used when --enable_readbase is enabled.'
            }
            if (!valid_postprocess_modes.contains(params.readbased_postprocess)) {
                exit 1, "ERROR: Unsupported --readbased_postprocess value '${params.readbased_postprocess}'. Valid options: all, phyloseq, taxburst, rgi_bwt."
            }
            if (params.skip_yacht) {
                exit 1, 'ERROR: --readbased_postprocess requires the merged Sourmash-YACHT result; disable --skip_yacht.'
            }
            if (params.readbased_postprocess == 'rgi_bwt' && !params.enable_rgi_bwt) {
                exit 1, 'ERROR: --enable_rgi_bwt is required for --readbased_postprocess rgi_bwt.'
            }
            if (params.readbased_postprocess == 'all' && !params.enable_rgi_bwt) {
                log.warn 'readbased_postprocess=all was requested, but RGI-BWT is disabled; the RGI-BWT post-processing output will be skipped while phyloseq and taxburst continue.'
            }
        }

        // Preflight: validate assembly_input usage
        if (params.assembly_input) {
            if (params.input_format != 'csv') {
                exit 1, "ERROR: --assembly_input requires --input_format to be 'csv'. Current value: ${params.input_format}"
            }
            if (params.enable_readbase) {
                exit 1, "ERROR: --assembly_input is not compatible with read-based profiling (--enable_readbase). Use assembly-based workflows only."
            }
            log.warn "⚠️  Pre-computed assemblies provided. Skipping preprocessing and assembly steps, jumping straight to binning."
            log.warn "⚠️  All preprocessing steps (fastp, hostile) will be skipped. Ensure input reads are already cleaned."
        }

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

        // Create input channel for reads
        def directory_pair_patterns = [
            "${params.input}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
            "${params.input}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}"
        ]
        input_ch = workflow.preview && !params.input ? channel.empty() : params.input_format == 'csv' ?
            createCsvInputChannel(params.input) :
            channel.fromFilePairs(directory_pair_patterns)
                .map { sample_id, reads ->
                    def meta = [id:sample_id, single_end:false, run_id:'default_run', group:'default_group']
                    return tuple(meta, reads)
                }

        // Handle pre-computed assemblies or preprocessing
        def cleaned_reads_source
        def assembly_megahit_contigs_ch = channel.empty()
        def assembly_metaspades_contigs_ch = channel.empty()
        def assembly_versions_ch = channel.empty()

        if (params.assembly_input) {
            // Use pre-computed assemblies - skip preprocessing
            cleaned_reads_source = input_ch
            
            // Parse assembly input CSV
            assembly_input_ch = createAssemblyInputChannel(params.assembly_input)
            
            // Separate compressed and uncompressed assemblies
            assembly_splits = assembly_input_ch
                .branch { _meta, assembly ->
                    compressed: assembly.name.endsWith('.gz')
                    uncompressed: !assembly.name.endsWith('.gz')
                }
            
            // Decompress gzipped assemblies
            GUNZIP_ASSEMBLIES(assembly_splits.compressed)
            assembly_decompressed = GUNZIP_ASSEMBLIES.out.gunzip
            
            // Combine decompressed and already uncompressed assemblies
            assembly_final_ch = assembly_decompressed.mix(assembly_splits.uncompressed)
            
            // All pre-computed assemblies go into megahit channel (used as generic assembly channel)
            // The actual assembler name is preserved in meta.assembler as metadata
            assembly_megahit_contigs_ch = assembly_final_ch
            
        } else {
            // Normal preprocessing path
            if (!params.skip_preprocess) {
                PREPROCESS(input_ch)
                cleaned_reads_source = PREPROCESS.out.cleaned_reads
            } else {
                // When skip_preprocess is true, use raw input reads directly
                cleaned_reads_source = input_ch
            }
        }

        // Initialize empty channels with default values
        def read_based_versions_ch = channel.empty()
        def read_based_results_ch = channel.empty()
        def read_based_rgi_ch = channel.empty()
        def read_based_postprocess_phyloseq_ch = channel.empty()
        def read_based_postprocess_taxburst_ch = channel.empty()
        def read_based_postprocess_rgi_bwt_ch = channel.empty()
        def binning_bam_ch = channel.empty()
        def binning_versions_ch = channel.empty()
        def binning_bins_ch = channel.empty()
        def binqc_versions_ch = channel.empty()
        def binqc_summary_ch = channel.empty()
        def binqc_quast_summary_ch = channel.empty()
        def binqc_multiqc_ch = channel.empty()
        def binannotation_versions_ch = channel.empty()
        def binannotation_multiqc_ch = channel.empty()
        def binclassification_versions_ch = channel.empty()
        def binclassification_summary_ch = channel.empty()

        // Use cleaned_reads_source for downstream workflows
        if (params.enable_readbase) {
            // Check if single-sample processing mode is enabled
            if (params.enable_singlesketch) {
                // Use READ_BASED_SINGLERUN for per-sample processing
                READ_BASED_SINGLERUN(cleaned_reads_source)
                read_based_versions_ch = READ_BASED_SINGLERUN.out.versions ?: channel.empty()
                read_based_results_ch = READ_BASED_SINGLERUN.out.results ?: channel.empty()
                read_based_rgi_ch = READ_BASED_SINGLERUN.out.rgi_results ?: channel.empty()
                read_based_postprocess_phyloseq_ch = READ_BASED_SINGLERUN.out.postprocess_phyloseq ?: channel.empty()
                read_based_postprocess_taxburst_ch = READ_BASED_SINGLERUN.out.postprocess_taxburst ?: channel.empty()
                read_based_postprocess_rgi_bwt_ch = READ_BASED_SINGLERUN.out.postprocess_rgi_bwt ?: channel.empty()
            } else {
                // Use READ_BASED for batch processing (default)
                READ_BASED(cleaned_reads_source)
                read_based_versions_ch = READ_BASED.out.versions ?: channel.empty()
                read_based_results_ch = READ_BASED.out.results ?: channel.empty()
                read_based_rgi_ch = READ_BASED.out.rgi_results ?: channel.empty()
                read_based_postprocess_phyloseq_ch = READ_BASED.out.postprocess_phyloseq ?: channel.empty()
                read_based_postprocess_taxburst_ch = READ_BASED.out.postprocess_taxburst ?: channel.empty()
                read_based_postprocess_rgi_bwt_ch = READ_BASED.out.postprocess_rgi_bwt ?: channel.empty()
            }
        } else {
            // Assembly-based pathway
            if (!params.assembly_input) {
                ASSEMBLY_BASED(cleaned_reads_source)
                assembly_versions_ch = ASSEMBLY_BASED.out.versions ?: channel.empty()
                assembly_megahit_contigs_ch = ASSEMBLY_BASED.out.megahit_contigs ?: channel.empty()
                assembly_metaspades_contigs_ch = ASSEMBLY_BASED.out.metaspades_contigs ?: channel.empty()
            }

            def all_assemblies = assembly_megahit_contigs_ch.mix(assembly_metaspades_contigs_ch)
            
            if (!params.skip_binning_bamabund) {
                BINNING_BAMABUND(all_assemblies, cleaned_reads_source)
                binning_bam_ch = BINNING_BAMABUND.out.bam_bai ?: channel.empty()
                binning_versions_ch = binning_versions_ch.mix(BINNING_BAMABUND.out.versions ?: channel.empty())

                // Run BINNING with outputs from BINNING_BAMABUND
                if (!params.skip_binning && !params.skip_binning_bamabund) {
                    BINNING(all_assemblies, BINNING_BAMABUND.out.bam_bai)
                    binning_versions_ch = binning_versions_ch.mix(BINNING.out.versions ?: channel.empty())

                    // Initialize channels for raw bins with metadata
                    def ch_raw_bins = BINNING.out.all_bins.map { meta, bin ->
                        def meta_new = meta.clone()
                        meta_new.refinement = 'unrefined'
                        [meta_new, bin]
                    }

                    // ============================================
                    // FIX: Run DOMAIN_CLASSIFICATION BEFORE refinement
                    // ============================================
                    if (params.bin_domain_classification && (params.skip_binning_refinement || params.postbinning_input != 'refined_bins_only')) {
                        // Group raw bins by assembly for domain classification
                        ch_bins_for_classification = ch_raw_bins
                            .map { meta, bin ->
                                def meta_clean = meta.clone()
                                meta_clean.remove('bin_id')
                                def group_key = "${meta_clean.id}_${meta_clean.assembler}"
                                [group_key, meta_clean, bin]
                            }
                            .groupTuple(by: 0)
                            .map { _key, metas, bins ->
                                [metas[0], bins.flatten()]
                            }

                        DOMAIN_CLASSIFICATION_RAW(all_assemblies, ch_bins_for_classification)
                        binning_versions_ch = binning_versions_ch.mix(DOMAIN_CLASSIFICATION_RAW.out.versions)

                        // Add domain metadata to raw bins
                        ch_domain_map = DOMAIN_CLASSIFICATION_RAW.out.classified_bins
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
                            .map { _bin_basename, meta, bin, domain ->
                                def meta_new = meta.clone()
                                meta_new.domain = domain ?: 'unclassified'
                                [meta_new, bin]
                            }
                    }

                    // ============================================
                    // Setup CheckM2 database (if needed)
                    // ============================================
                    // CheckM2 database is only required if:
                    // - Using CheckM2 for bin QC (binqc_tool == 'checkm2'), OR
                    // - Using Binette for bin refinement (refine_tool == 'binette')
                    //
                    // Other combinations like CheckM + DAS_TOOL do NOT need CheckM2 database
                    // (e.g., CheckM = legacy checkm v1, DAS_TOOL = only needs contig2bin tables)
                    //
                    def ch_checkm2_db = channel.empty()
                    def needs_checkm2_db = (params.binqc_tool == 'checkm2' && !params.skip_binqc) || 
                                           (params.refine_tool == 'binette' && !params.skip_binning_refinement)
                    
                    if (needs_checkm2_db && !params.checkm2_db) {
                        // Download CheckM2 database automatically if needed and not provided
                        log.info "Downloading CheckM2 database for ${params.binqc_tool == 'checkm2' ? 'CheckM2 QC' : ''}${(params.binqc_tool == 'checkm2' && params.refine_tool == 'binette') ? ' and ' : ''}${params.refine_tool == 'binette' ? 'Binette refinement' : ''}"
                        CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_version)
                        // Extract just the database path from the tuple (meta, path) and
                        // convert it to a value channel (single shared DB for all downstream tasks)
                        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
                                                        .map { _meta, db -> db }
                        binning_versions_ch = binning_versions_ch.mix(CHECKM2_DATABASEDOWNLOAD.out.versions_checkm2_databasedownload)
                    } else if (params.checkm2_db) {
                        // Use provided CheckM2 database - pass just the path
                        ch_checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
                    } else {
                        // No CheckM2 database needed - create empty channel with null value
                        // This safely handles cases like: CheckM (v1) + DAS_TOOL, BUSCO + DAS_TOOL, etc.
                        ch_checkm2_db = channel.value([])
                    }

                    // ============================================
                    // Run BINNING_REFINEMENT (with domain already set)
                    // ============================================
                    def ch_refined_bins = channel.empty()

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

                        BINNING_REFINEMENT(all_assemblies, ch_bins_for_refinement, ch_checkm2_db)
                        
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

                        // Classify refined bins separately because refined bin names
                        // do not match the original raw-bin basenames.
                        if (params.bin_domain_classification && params.postbinning_input != 'raw_bins_only') {
                            ch_refined_bins_for_classification = ch_refined_bins
                                .map { meta, bin ->
                                    def meta_clean = meta.clone()
                                    meta_clean.remove('bin_id')
                                    def group_key = "${meta_clean.id}_${meta_clean.assembler}"
                                    [group_key, meta_clean, bin]
                                }
                                .groupTuple(by: 0)
                                .map { _key, metas, bins ->
                                    [metas[0], bins.flatten()]
                                }

                            DOMAIN_CLASSIFICATION_REFINED(all_assemblies, ch_refined_bins_for_classification)
                            ch_refined_bins = DOMAIN_CLASSIFICATION_REFINED.out.classified_bins
                            binning_versions_ch = binning_versions_ch.mix(DOMAIN_CLASSIFICATION_REFINED.out.versions)
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

                        BIN_QC(ch_bins_for_qc, ch_checkm2_db)
                        binqc_versions_ch = BIN_QC.out.versions ?: channel.empty()
                        binqc_summary_ch = BIN_QC.out.qc_summary ?: channel.empty()
                        binqc_quast_summary_ch = BIN_QC.out.quast_summary ?: channel.empty()
                        binqc_multiqc_ch = BIN_QC.out.multiqc_files ?: channel.empty()
                        binning_versions_ch = binning_versions_ch.mix(binqc_versions_ch)
                    }

                    // Run BIN_ANNOTATION in parallel with BIN_QC
                    if (!params.skip_bin_annotation && !params.skip_binning) {
                        ch_bins_for_annotation = ch_bins_split.for_annotation

                        BIN_ANNOTATION(ch_bins_for_annotation)
                        binannotation_versions_ch = BIN_ANNOTATION.out.versions ?: channel.empty()
                        binannotation_multiqc_ch = BIN_ANNOTATION.out.multiqc_files ?: channel.empty()
                        binning_versions_ch = binning_versions_ch.mix(binannotation_versions_ch)
                    }

                    // Run BIN_CLASSIFICATION in parallel with BIN_QC and BIN_ANNOTATION
                    if (!params.skip_bin_classification && !params.skip_binning) {
                        ch_bins_for_classification = ch_bins_split.for_classification

                        BIN_CLASSIFICATION(ch_bins_for_classification)
                        binclassification_versions_ch = BIN_CLASSIFICATION.out.versions ?: channel.empty()
                        binclassification_summary_ch = BIN_CLASSIFICATION.out.genome_classification ?: channel.empty()
                        binning_versions_ch = binning_versions_ch.mix(binclassification_versions_ch)
                    }
                }
            }
        }

        def combined_versions = channel.empty()
        combined_versions = combined_versions.mix(assembly_versions_ch ?: channel.empty())
        combined_versions = combined_versions.mix(read_based_versions_ch ?: channel.empty())
        combined_versions = combined_versions.mix(binning_versions_ch ?: channel.empty())

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
            cleanup_sources = cleanup_sources.findAll { source -> source }
            if (cleanup_sources) {
                def cleanup_trigger = cleanup_sources[0]
                if (cleanup_sources.size() > 1) {
                    cleanup_sources.drop(1).each { ch ->
                        cleanup_trigger = cleanup_trigger.mix(ch)
                    }
                }
                CLEANUP(cleanup_trigger.collect())
            }
        }

    emit:
        versions = combined_versions.ifEmpty(null)
        results = read_based_results_ch.ifEmpty(null)
        read_based_rgi = read_based_rgi_ch.ifEmpty(null)
        read_based_postprocess_phyloseq = read_based_postprocess_phyloseq_ch.ifEmpty(null)
        read_based_postprocess_taxburst = read_based_postprocess_taxburst_ch.ifEmpty(null)
        read_based_postprocess_rgi_bwt = read_based_postprocess_rgi_bwt_ch.ifEmpty(null)
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

workflow {
    if (!(workflow.preview && !params.input && !params.assembly_input)) {
        METAFLOW()
    }
}
