#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MEGAHIT } from '../../modules/nf-core/megahit/main'
include { SPADES as METASPADES } from '../../modules/nf-core/spades/main'
include { QUAST_METAQUAST as METAQUAST } from '../../modules/local/quast/metaquast/main'

workflow ASSEMBLY_BASED {
    take:
    cleaned_reads_ch // Channel: [ val(meta), [ reads ] ]

    main:
    versions_ch = Channel.empty()
    ch_megahit_outputs = Channel.empty()
    ch_metaspades_outputs = Channel.empty()
    ch_assemblies_for_quast = Channel.empty()

    // MEGAHIT Assembly
    if (!params.skip_megahit) {
        ch_megahit_input = cleaned_reads_ch.map { meta, reads -> 
            def reads_list = reads.collect()
            def sorted_reads = reads_list.sort()
            def read1 = sorted_reads[0]
            def read2 = sorted_reads.size() > 1 ? sorted_reads[1] : null
            [meta, read1 ? [read1] : [], read2 ? [read2] : []]
        }

        MEGAHIT(ch_megahit_input)
        versions_ch = versions_ch.mix(MEGAHIT.out.versions)
        ch_megahit_outputs = MEGAHIT.out.contigs.map { meta, assembly ->
            def meta_new = meta + [assembler: 'MEGAHIT']
            [meta_new, assembly]
        }
        ch_assemblies_for_quast = ch_assemblies_for_quast.mix(ch_megahit_outputs)
    }

    // METASPADES Assembly
    if (!params.skip_spades) {
        ch_spades_input = cleaned_reads_ch.map { meta, reads -> 
            def reads_list = reads.collect()
            def sorted_reads = reads_list.sort()
            def illumina = sorted_reads.size() > 1 ? sorted_reads : [sorted_reads[0]]
            meta.single_end = sorted_reads.size() <= 1
            [meta, illumina, [], []]
        }

        METASPADES(ch_spades_input, [], [])
        versions_ch = versions_ch.mix(METASPADES.out.versions)
        ch_metaspades_outputs = METASPADES.out.contigs.map { meta, assembly ->
            def meta_new = meta + [assembler: 'SPAdes']
            [meta_new, assembly]
        }
        ch_assemblies_for_quast = ch_assemblies_for_quast.mix(ch_metaspades_outputs)
    }

    // QUAST Analysis
    if (!params.skip_quast) {
        METAQUAST(ch_assemblies_for_quast)
        versions_ch = versions_ch.mix(METAQUAST.out.versions)
    }

    emit:
    megahit_contigs = ch_megahit_outputs
    metaspades_contigs = ch_metaspades_outputs
    versions = versions_ch.ifEmpty(null)
    quast_qc = params.skip_quast ? null : METAQUAST.out.qc
}
