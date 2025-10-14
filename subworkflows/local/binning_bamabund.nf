#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Include nf-core modules
include { BOWTIE2_BUILD as BOWTIE2_ASSEMBLY_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ASSEMBLY_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { BWAMEM2_INDEX as BWAMEM2_ASSEMBLY_INDEX } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_ASSEMBLY_MEM }     from '../../modules/nf-core/bwamem2/mem/main'
// include { SAMTOOLS_VIEW }                           from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT }                           from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX }                          from '../../modules/nf-core/samtools/index/main'

workflow BINNING_BAMABUND {
    take:
    assemblies      // channel: [ val(meta), path(contigs) ] - from MEGAHIT/METASPADES
    cleaned_reads   // channel: [ val(meta), [ reads ] ] - from HOSTILE_CLEAN

    main:
    versions_ch = Channel.empty()
    
    // Prepare assemblies channel by decompressing if needed
    ch_assemblies_for_mapping = assemblies.map { meta, assembly ->
        def new_meta = meta.clone()
        new_meta.assembler = meta.assembler ?: 'unknown'
        [new_meta, assembly]
    }

    // Match assemblies with their corresponding reads
    // Group by sample ID to ensure correct pairing
    ch_assemblies_by_sample = ch_assemblies_for_mapping.map { meta, assembly ->
        [meta.id, meta, assembly]
    }
    
    ch_reads_by_sample = cleaned_reads.map { meta, reads ->
        [meta.id, meta, reads]
    }
    
    ch_assembly_reads_paired = ch_assemblies_by_sample
        .combine(ch_reads_by_sample, by: 0)
        .map { sample_id, assembly_meta, assembly, reads_meta, reads ->
            // Merge metadata, prioritizing assembly metadata
            def merged_meta = assembly_meta.clone()
            [merged_meta, assembly, reads]
        }

    // Branch based on mapper tool selection
    if (params.binprepare_mappertool == 'bwamem2') {
        // BWAMEM2 workflow
        BWAMEM2_ASSEMBLY_INDEX(
            ch_assembly_reads_paired.map { meta, assembly, reads -> [meta, assembly] }
        )
        versions_ch = versions_ch.mix(BWAMEM2_ASSEMBLY_INDEX.out.versions)

        // Prepare input for BWAMEM2_MEM: [meta, reads], [meta, index], [meta, fasta]
        ch_bwamem2_input = ch_assembly_reads_paired
            .map { meta, assembly, reads -> [meta.id, meta, reads, assembly] }
            .combine(BWAMEM2_ASSEMBLY_INDEX.out.index.map { meta, index -> [meta.id, meta, index] }, by: 0)
            .map { sample_id, reads_meta, reads, assembly, index_meta, index ->
                [reads_meta, reads, index, assembly]
            }

        BWAMEM2_ASSEMBLY_MEM(
            ch_bwamem2_input.map { meta, reads, index, fasta -> [meta, reads] },
            ch_bwamem2_input.map { meta, reads, index, fasta -> [meta, index] },
            ch_bwamem2_input.map { meta, reads, index, fasta -> [meta, fasta] },
            false  // sort_bam = false, we'll sort separately
        )
        versions_ch = versions_ch.mix(BWAMEM2_ASSEMBLY_MEM.out.versions)

        ch_mapped_bam = BWAMEM2_ASSEMBLY_MEM.out.bam
        
    } else {
        // BOWTIE2 workflow (default)
        BOWTIE2_ASSEMBLY_BUILD(
            ch_assembly_reads_paired.map { meta, assembly, reads -> [meta, assembly] }
        )
        versions_ch = versions_ch.mix(BOWTIE2_ASSEMBLY_BUILD.out.versions)

        // Prepare input for BOWTIE2_ALIGN: [meta, reads], [meta, index], [meta, fasta]
        ch_bowtie2_input = ch_assembly_reads_paired
            .map { meta, assembly, reads -> [meta.id, meta, reads, assembly] }
            .combine(BOWTIE2_ASSEMBLY_BUILD.out.index.map { meta, index -> [meta.id, meta, index] }, by: 0)
            .map { sample_id, reads_meta, reads, assembly, index_meta, index ->
                [reads_meta, reads, index, assembly]
            }

        BOWTIE2_ASSEMBLY_ALIGN(
            ch_bowtie2_input.map { meta, reads, index, fasta -> [meta, reads] },
            ch_bowtie2_input.map { meta, reads, index, fasta -> [meta, index] },
            ch_bowtie2_input.map { meta, reads, index, fasta -> [meta, fasta] },
            false,  // save_unaligned = false
            false   // sort_bam = false, we'll sort separately
        )
        versions_ch = versions_ch.mix(BOWTIE2_ASSEMBLY_ALIGN.out.versions)

        ch_mapped_bam = BOWTIE2_ASSEMBLY_ALIGN.out.bam
    }

    // Sort BAM files
    SAMTOOLS_SORT(
        ch_mapped_bam,
        [[id: 'null'], []],  // No reference fasta needed for sorting
        'bai'  // Create BAI index
    )
    versions_ch = versions_ch.mix(SAMTOOLS_SORT.out.versions)

    // Index sorted BAM files (if not already indexed by SAMTOOLS_SORT)
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    versions_ch = versions_ch.mix(SAMTOOLS_INDEX.out.versions)

    // Prepare output channels for downstream binning
    // Combine BAM with BAI index
    ch_sorted_bam_bai = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: 0)
        .map { meta, bam, bai ->
            [meta, bam, bai]
        }

    emit:
    bam         = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), path(bam) ]
    bai         = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), path(bai) ]
    bam_bai     = ch_sorted_bam_bai               // channel: [ val(meta), path(bam), path(bai) ]
    versions    = versions_ch
}
