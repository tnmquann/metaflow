#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Include nf-core modules
include { BOWTIE2_BUILD as BOWTIE2_ASSEMBLY_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ASSEMBLY_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { BWAMEM2_INDEX as BWAMEM2_ASSEMBLY_INDEX } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_ASSEMBLY_MEM }     from '../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT }                           from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX }                          from '../../modules/nf-core/samtools/index/main'

workflow BINNING_BAMABUND {
    take:
    assemblies      // channel: [ val(meta), path(contigs) ]
    cleaned_reads   // channel: [ val(meta), [ reads ] ]

    main:
    versions_ch = Channel.empty()
    
    // Prepare assemblies channel
    ch_assemblies_for_mapping = assemblies.map { meta, assembly ->
        def new_meta = meta.clone()
        new_meta.assembler = meta.assembler ?: 'unknown'
        [new_meta, assembly]
    }

    // Match assemblies with their corresponding reads
    ch_assemblies_by_sample = ch_assemblies_for_mapping.map { meta, assembly ->
        [meta.id, meta, assembly]
    }
    
    ch_reads_by_sample = cleaned_reads.map { meta, reads ->
        [meta.id, meta, reads]
    }
    
    ch_assembly_reads_paired = ch_assemblies_by_sample
        .combine(ch_reads_by_sample, by: 0)
        .map { sample_id, assembly_meta, assembly, reads_meta, reads ->
            // Create a unique identifier that includes BOTH assembly and reads info
            def merged_meta = assembly_meta.clone()
            merged_meta.reads_id = reads_meta.id  // Add reads ID to distinguish alignments
            [merged_meta, assembly, reads]
        }

    // Branch based on mapper tool selection
    if (params.binprepare_mappertool == 'bwamem2') {
        // BWAMEM2 workflow
        BWAMEM2_ASSEMBLY_INDEX(
            ch_assembly_reads_paired.map { meta, assembly, reads -> [meta, assembly] }
        )
        versions_ch = versions_ch.mix(BWAMEM2_ASSEMBLY_INDEX.out.versions)

        ch_bwamem2_input = ch_assembly_reads_paired
            .map { meta, assembly, reads -> [meta.id + "_" + meta.assembler, meta, reads, assembly] }
            .combine(BWAMEM2_ASSEMBLY_INDEX.out.index.map { meta, index -> [meta.id + "_" + meta.assembler, meta, index] }, by: 0)
            .map { key, reads_meta, reads, assembly, index_meta, index ->
                [reads_meta, reads, index, assembly]
            }

        BWAMEM2_ASSEMBLY_MEM(
            ch_bwamem2_input.map { meta, reads, index, fasta -> [meta, reads] },
            ch_bwamem2_input.map { meta, reads, index, fasta -> [meta, index] },
            ch_bwamem2_input.map { meta, reads, index, fasta -> [meta, fasta] },
            false
        )
        versions_ch = versions_ch.mix(BWAMEM2_ASSEMBLY_MEM.out.versions)

        ch_mapped_bam = BWAMEM2_ASSEMBLY_MEM.out.bam
        
    } else {
        // BOWTIE2 workflow (default)
        BOWTIE2_ASSEMBLY_BUILD(
            ch_assembly_reads_paired
                .map { meta, assembly, reads -> [meta, assembly] }
                .unique { it[0].id + "_" + it[0].assembler }  // Only build index once per assembly
        )
        versions_ch = versions_ch.mix(BOWTIE2_ASSEMBLY_BUILD.out.versions)

        ch_bowtie2_input = ch_assembly_reads_paired
            .map { meta, assembly, reads -> [meta.id + "_" + meta.assembler, meta, reads, assembly] }
            .combine(BOWTIE2_ASSEMBLY_BUILD.out.index.map { meta, index -> [meta.id + "_" + meta.assembler, meta, index] }, by: 0)
            .map { key, reads_meta, reads, assembly, index_meta, index ->
                [reads_meta, reads, index, assembly]
            }

        BOWTIE2_ASSEMBLY_ALIGN(
            ch_bowtie2_input.map { meta, reads, index, fasta -> [meta, reads] },
            ch_bowtie2_input.map { meta, reads, index, fasta -> [meta, index] },
            ch_bowtie2_input.map { meta, reads, index, fasta -> [meta, fasta] },
            false,
            false
        )
        versions_ch = versions_ch.mix(BOWTIE2_ASSEMBLY_ALIGN.out.versions)

        ch_mapped_bam = BOWTIE2_ASSEMBLY_ALIGN.out.bam
    }

    // Sort BAM files
    SAMTOOLS_SORT(
        ch_mapped_bam,
        [[id: 'null'], []],
        'bai'
    )
    versions_ch = versions_ch.mix(SAMTOOLS_SORT.out.versions)

    // Index sorted BAM files
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    versions_ch = versions_ch.mix(SAMTOOLS_INDEX.out.versions)

    // Group BAM/BAI by assembly (removing reads_id from grouping key)
    ch_sorted_bam_bai = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: 0)
        .map { meta, bam, bai ->
            // Create grouping key without reads_id
            def group_meta = meta.clone()
            group_meta.remove('reads_id')
            [group_meta.id + "_" + group_meta.assembler, group_meta, bam, bai]
        }
        .groupTuple(by: 0)
        .map { key, metas, bams, bais ->
            // Use the first meta as representative
            [metas[0], bams, bais]
        }

    emit:
    bam         = SAMTOOLS_SORT.out.bam
    bai         = SAMTOOLS_INDEX.out.bai
    bam_bai     = ch_sorted_bam_bai
    versions    = versions_ch
}
