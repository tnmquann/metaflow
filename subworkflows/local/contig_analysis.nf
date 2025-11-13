#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import modules
include { SKANI_SEARCH } from '../../modules/nf-core/skani/search/main'
include { TIARA_TIARA } from '../../modules/nf-core/tiara/tiara/main'
include { SEQKIT_FX2TAB } from '../../modules/nf-core/seqkit/fx2tab/main'
include { GENOMAD_ENDTOEND } from '../../modules/nf-core/genomad/endtoend/main'
include { GENOMAD_DOWNLOAD } from '../../modules/nf-core/genomad/download/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_CONTIGS } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_CONTIGS } from '../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_COVERAGE } from '../../modules/local/samtools/coverage/main'

// Annotation modules
include { PROKKA } from '../../modules/nf-core/prokka/main'
include { PYRODIGAL } from '../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FNA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FAA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_GBK } from '../../modules/nf-core/gunzip/main'

workflow CONTIG_ANALYSIS {
    take:
    ch_assemblies      // Channel: [ val(meta), path(contigs) ] - from MEGAHIT/METASPADES
    ch_cleaned_reads   // Channel: [ val(meta), [ reads ] ] - from HOSTILE_CLEAN

    main:
    ch_versions = Channel.empty()

    // Step 1: Index assemblies with BWA-MEM2
    BWAMEM2_INDEX_CONTIGS(ch_assemblies)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX_CONTIGS.out.versions)

    // Step 2: Prepare input for BWA-MEM2 MEM (join reads with index by meta.id)
    ch_reads_with_index = ch_cleaned_reads
        .map { meta, reads -> 
            def key = meta.id
            [key, meta, reads]
        }
        .combine(
            BWAMEM2_INDEX_CONTIGS.out.index.map { meta, index ->
                def key = meta.id
                [key, meta, index]
            },
            by: 0
        )
        .map { key, meta_reads, reads, meta_assembly, index ->
            // Merge meta information, keeping assembler from assembly
            def meta_merged = meta_reads + [assembler: meta_assembly.assembler]
            [meta_merged, reads, meta_assembly, index]
        }

    // Step 3: Map reads to assemblies (sort BAM output)
    BWAMEM2_MEM_CONTIGS(
        ch_reads_with_index.map { meta, reads, meta_asm, index -> [meta, reads] },
        ch_reads_with_index.map { meta, reads, meta_asm, index -> [meta_asm, index] },
        ch_assemblies.map { meta, assembly -> [meta, assembly] },
        true  // sort_bam = true
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM_CONTIGS.out.versions)

    // Step 4: Calculate coverage statistics (simplified - only BAM input)
    SAMTOOLS_COVERAGE(BWAMEM2_MEM_CONTIGS.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    // Step 5: Extract contig sequence metadata
    SEQKIT_FX2TAB(ch_assemblies)
    ch_versions = ch_versions.mix(SEQKIT_FX2TAB.out.versions)

    // Step 6 (optional): TIARA domain classification
    ch_tiara_results = Channel.empty()
    if (!params.skip_tiara_contigqc) {
        TIARA_TIARA(ch_assemblies)
        ch_tiara_results = TIARA_TIARA.out.classifications
        ch_versions = ch_versions.mix(TIARA_TIARA.out.versions)
    }

    // Step 7 (optional): SKANI search for ANI-based classification
    ch_skani_results = Channel.empty()
    if (!params.skip_skani_search) {
        if (params.skani_db) {
            ch_skani_db = Channel.value([
                [id: 'skani_db'], 
                file(params.skani_db, checkIfExists: true)
            ])
            SKANI_SEARCH(ch_assemblies, ch_skani_db)
            ch_skani_results = SKANI_SEARCH.out.search
            ch_versions = ch_versions.mix(SKANI_SEARCH.out.versions)
        } else {
            log.warn "SKANI_SEARCH: No database provided (--skani_db). Skipping SKANI analysis."
        }
    }

    // Step 8 (optional): geNomad for plasmid/virus identification
    ch_genomad_plasmid_summary = Channel.empty()
    ch_genomad_virus_summary = Channel.empty()
    if (!params.skip_genomad_e2e) {
        // Prepare or download geNomad database
        if (params.genomad_db) {
            ch_genomad_db = file(params.genomad_db, checkIfExists: true)
        } else {
            GENOMAD_DOWNLOAD()
            ch_genomad_db = GENOMAD_DOWNLOAD.out.genomad_db
            ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions)
            
            if (params.save_genomad_db) {
                log.info "geNomad database downloaded and will be saved to output directory."
            }
        }

        GENOMAD_ENDTOEND(ch_assemblies, ch_genomad_db)
        ch_genomad_plasmid_summary = GENOMAD_ENDTOEND.out.plasmid_summary
        ch_genomad_virus_summary = GENOMAD_ENDTOEND.out.virus_summary
        ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions)
    }

    // Step 9 (optional): Contig Annotation
    ch_annotation_faa = Channel.empty()
    ch_annotation_fna = Channel.empty()
    ch_annotation_gbk = Channel.empty()
    
    if (!params.skip_contig_annotation) {
        if (params.annotation_tool == 'pyrodigal') {
            PYRODIGAL(ch_assemblies, "gbk")
            
            GUNZIP_PYRODIGAL_FAA(PYRODIGAL.out.faa)
            GUNZIP_PYRODIGAL_FNA(PYRODIGAL.out.fna)
            GUNZIP_PYRODIGAL_GBK(PYRODIGAL.out.annotations)
            
            ch_annotation_faa = GUNZIP_PYRODIGAL_FAA.out.gunzip
            ch_annotation_fna = GUNZIP_PYRODIGAL_FNA.out.gunzip
            ch_annotation_gbk = GUNZIP_PYRODIGAL_GBK.out.gunzip
            
            ch_versions = ch_versions.mix(PYRODIGAL.out.versions)
            ch_versions = ch_versions.mix(GUNZIP_PYRODIGAL_FAA.out.versions)
            ch_versions = ch_versions.mix(GUNZIP_PYRODIGAL_FNA.out.versions)
            ch_versions = ch_versions.mix(GUNZIP_PYRODIGAL_GBK.out.versions)
            
        } else if (params.annotation_tool == 'prokka') {
            PROKKA(ch_assemblies, [], [])
            
            ch_annotation_faa = PROKKA.out.faa
            ch_annotation_fna = PROKKA.out.fna
            ch_annotation_gbk = PROKKA.out.gbk
            
            ch_versions = ch_versions.mix(PROKKA.out.versions)
        } else {
            log.warn "Invalid annotation_tool: ${params.annotation_tool}. Skipping annotation. Valid options: 'pyrodigal', 'prokka'"
        }
    }

    // Step 10: Combine QC results for downstream analysis
    ch_contig_qc_bundle = SAMTOOLS_COVERAGE.out.coverage
        .join(SEQKIT_FX2TAB.out.text, by: 0, remainder: true)

    // Add optional outputs if they exist
    if (!params.skip_skani_search && params.skani_db) {
        ch_contig_qc_bundle = ch_contig_qc_bundle
            .join(ch_skani_results, by: 0, remainder: true)
    }

    if (!params.skip_genomad_e2e) {
        ch_contig_qc_bundle = ch_contig_qc_bundle
            .join(ch_genomad_plasmid_summary, by: 0, remainder: true)
    }

    if (!params.skip_tiara_contigqc) {
        ch_contig_qc_bundle = ch_contig_qc_bundle
            .join(ch_tiara_results, by: 0, remainder: true)
    }

    emit:
    bam                = BWAMEM2_MEM_CONTIGS.out.bam
    coverage           = SAMTOOLS_COVERAGE.out.coverage
    seqkit_stats       = SEQKIT_FX2TAB.out.text
    tiara_class        = ch_tiara_results
    skani_results      = ch_skani_results
    genomad_plasmids   = ch_genomad_plasmid_summary
    genomad_viruses    = ch_genomad_virus_summary
    annotation_faa     = ch_annotation_faa
    annotation_fna     = ch_annotation_fna
    annotation_gbk     = ch_annotation_gbk
    contig_qc_bundle   = ch_contig_qc_bundle
    versions           = ch_versions
}