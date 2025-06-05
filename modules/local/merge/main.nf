process MERGE_PAIREDENDSEQS {
    tag { meta.id }
    label 'process_low'
    maxForks 2  // Allow 2 parallel executions

    conda "${projectDir}/env/read_based.yaml"

    publishDir "${params.merged_seq_dir}", mode: params.publish_dir_mode, enabled: params.enable_copymergedseqs, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/${meta.id}.fastq.gz"), emit: merged_seqs
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    mkdir -p ${meta.id}
    cat ${reads} > ${meta.id}/${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 || echo "cat (coreutils)")
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${prefix}.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: "stub_version"
    END_VERSIONS
    """
}