process MERGE_SEQUENCES {
    tag "Merging sequences"
    publishDir "${params.base_dir}/merged_seq",
        mode: 'copy',
        enabled: params.enable_copymergedseqs

    input:
    path trimmed_fastq

    output:
    path "merged_seq", emit: merged_seqs

    script:
    """
    mkdir -p merged_seq
    for file in ${trimmed_fastq}/*_1.fastq.gz; do
        base=\$(basename \$file _1.fastq.gz)
        cat ${trimmed_fastq}/\${base}_1.fastq.gz ${trimmed_fastq}/\${base}_2.fastq.gz > merged_seq/\${base}.fastq.gz
    done
    """
}