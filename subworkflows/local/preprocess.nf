#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import nf-core modules
include { FASTQC as FASTQC_RAW     } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_RMHOST  } from '../../modules/nf-core/fastqc/main'
include { FASTP                     } from '../../modules/nf-core/fastp/main'
include { HOSTILE_FETCH             } from '../../modules/nf-core/hostile/fetch/main'
include { HOSTILE_CLEAN             } from '../../modules/nf-core/hostile/clean/main'

workflow PREPROCESS {
    take:
    reads_ch // channel: [ val(meta), path(reads) ]

    main:
    versions_coll = Channel.empty()

    // 1. Raw QC
    FASTQC_RAW ( reads_ch )
    versions_coll = versions_coll.mix(FASTQC_RAW.out.versions)

    // 2. Trimming with Fastp
    FASTP (
        reads_ch,
        params.fastp_adapter_fasta ?: [],   // Empty list if null
        false,                              // discard_trimmed_pass
        params.fastp_save_trimmed_fail ?: false,
        params.fastp_save_merged ?: false
    )
    versions_coll = versions_coll.mix(FASTP.out.versions)

    // 3. Trimmed QC
    FASTQC_TRIMMED ( FASTP.out.reads )
    versions_coll = versions_coll.mix(FASTQC_TRIMMED.out.versions)

    // 4. Host Removal
    ch_hostile_ref_dir = Channel.empty()
    if (params.hostile_reference) {
        ch_hostile_ref_dir = Channel.value(file(params.hostile_reference, checkIfExists: true))
    } else {
        // Use the hostile_index parameter as input for HOSTILE_FETCH
        def index_name = params.hostile_index ?: 'human-t2t-hla'
        HOSTILE_FETCH ( index_name )
        ch_hostile_ref_dir = HOSTILE_FETCH.out.reference
        versions_coll = versions_coll.mix(HOSTILE_FETCH.out.versions)
    }

    // Before HOSTILE_CLEAN process call, prepare the reference channel properly
    def ref_name = params.hostile_ref_name ?: params.hostile_index
    HOSTILE_CLEAN (
        FASTP.out.reads,
        ch_hostile_ref_dir.map { ref_dir -> [ ref_name, ref_dir ] }
    )
    versions_coll = versions_coll.mix(HOSTILE_CLEAN.out.versions)

    // 5. Host-Removed QC
    FASTQC_RMHOST ( HOSTILE_CLEAN.out.fastq )
    versions_coll = versions_coll.mix(FASTQC_RMHOST.out.versions)

    emit:
    cleaned_reads      = HOSTILE_CLEAN.out.fastq
    versions           = versions_coll.unique().collect()
    fastqc_raw_zip     = FASTQC_RAW.out.zip
    fastqc_trimmed_zip = FASTQC_TRIMMED.out.zip
    fastqc_rmhost_zip  = FASTQC_RMHOST.out.zip
    fastp_reports      = FASTP.out.json.join(FASTP.out.html).join(FASTP.out.log)
    hostile_clean_json = HOSTILE_CLEAN.out.json
}