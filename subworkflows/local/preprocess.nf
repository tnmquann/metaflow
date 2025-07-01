#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import nf-core modules from your ./modules/nf-core directory
include { FASTQC as FASTQC_RAW     } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_RMHOST  } from '../../modules/nf-core/fastqc/main'
// include { MULTIQC as MULTIQC_RAW     } from '../../modules/nf-core/multiqc/main'
// include { MULTIQC as MULTIQC_TRIMMED } from '../../modules/nf-core/multiqc/main'
// include { MULTIQC as MULTIQC_RMHOST  } from '../../modules/nf-core/multiqc/main'
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

    // MULTIQC_RAW (
    //     FASTQC_RAW.out.zip.collect(),
    //     params.multiqc_config ?: [],        // Empty list if null
    //     params.multiqc_extra_config ?: [],  // Empty list if null
    //     params.multiqc_logo ?: [],          // Empty list if null
    //     params.multiqc_replace_names ?: [], // Empty list if null
    //     params.multiqc_sample_names ?: []   // Empty list if null
    // )
    // versions_coll = versions_coll.mix(MULTIQC_RAW.out.versions)

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

    // MULTIQC_TRIMMED (
    //     FASTP.out.json.mix(FASTQC_TRIMMED.out.zip).collect(),
    //     params.multiqc_config ?: [],
    //     params.multiqc_extra_config ?: [],
    //     params.multiqc_logo ?: [],
    //     params.multiqc_replace_names ?: [],
    //     params.multiqc_sample_names ?: []
    // )
    // versions_coll = versions_coll.mix(MULTIQC_TRIMMED.out.versions)

    // 4. Host Removal
    ch_hostile_ref_dir = Channel.empty()
    if (params.hostile_reference) {
        ch_hostile_ref_dir = Channel.value(file(params.hostile_reference, checkIfExists: true))
    } else {
        HOSTILE_FETCH ( )
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

    // MULTIQC_RMHOST (
    //     HOSTILE_CLEAN.out.json.mix(FASTQC_RMHOST.out.zip).collect(),
    //     params.multiqc_config ?: [],
    //     params.multiqc_extra_config ?: [],
    //     params.multiqc_logo ?: [],
    //     params.multiqc_replace_names ?: [],
    //     params.multiqc_sample_names ?: []
    // )
    // versions_coll = versions_coll.mix(MULTIQC_RMHOST.out.versions)

    emit:
    cleaned_reads      = HOSTILE_CLEAN.out.fastq
    versions           = versions_coll.unique().collect()

    fastqc_raw_zip        = FASTQC_RAW.out.zip
    // multiqc_raw_report    = MULTIQC_RAW.out.report
    fastqc_trimmed_zip    = FASTQC_TRIMMED.out.zip
    // multiqc_trimmed_report= MULTIQC_TRIMMED.out.report
    fastqc_rmhost_zip     = FASTQC_RMHOST.out.zip
    // multiqc_rmhost_report = MULTIQC_RMHOST.out.report
    fastp_reports         = FASTP.out.json.join(FASTP.out.html).join(FASTP.out.log)
    hostile_clean_json    = HOSTILE_CLEAN.out.json
}