process MERGE_SEQUENCES {
    tag "Merging sequences"
    label 'process_low'

    conda "${projectDir}/env/read_based.yaml"

    publishDir "${params.merged_seq_dir}",
        mode: params.publish_dir_mode,
        enabled: params.enable_copymergedseqs,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    path trimmed_fastq    // Directory containing paired FASTQ files

    output:
    path "merged_seq", emit: merged_seqs
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p merged_seq
    for file in ${trimmed_fastq}/*_1.fastq.gz; do
        base=\$(basename \$file _1.fastq.gz)
        cat ${trimmed_fastq}/\${base}_1.fastq.gz ${trimmed_fastq}/\${base}_2.fastq.gz > merged_seq/\${base}.fastq.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 || echo "cat (coreutils)")
        basename: \$(basename --version 2>&1 | head -n 1 || echo "basename (coreutils)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p merged_seq
    touch merged_seq/test.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: "stub_version"
        basename: "stub_version"
    END_VERSIONS
    """
}