#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MEGAHIT } from '../../modules/nf-core/megahit/main'
include { SPADES as METASPADES } from '../../modules/nf-core/spades/main'

workflow ASSEMBLY_BASED {
    take:
    cleaned_reads_ch // Channel: [ val(meta), [ reads ] ]

    main:
    versions_ch = Channel.empty()
    ch_megahit_outputs = Channel.empty()
    ch_metaspades_outputs = Channel.empty()

    // MEGAHIT Assembly - Keep existing implementation
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
        ch_megahit_outputs = MEGAHIT.out.contigs
    }

    // METASPADES Assembly - Modified implementation  
    if (!params.skip_spades) {
        ch_spades_input = cleaned_reads_ch.map { meta, reads -> 
            def reads_list = reads.collect()
            def sorted_reads = reads_list.sort()
            def illumina = sorted_reads.size() > 1 ? sorted_reads : [sorted_reads[0]]
            // Ensure meta.single_end is properly set
            meta.single_end = sorted_reads.size() <= 1
            [meta, illumina, [], []] // illumina, pacbio, nanopore
        }

        METASPADES(
            ch_spades_input,
            [],  // Empty yml channel
            []   // Empty hmm channel
        )
        
        versions_ch = versions_ch.mix(METASPADES.out.versions)
        // Prefer contigs over scaffolds for downstream analysis
        ch_metaspades_outputs = METASPADES.out.contigs
    }

    emit:
    megahit_contigs = ch_megahit_outputs
    metaspades_contigs = ch_metaspades_outputs
    versions = versions_ch.ifEmpty(null)
}
