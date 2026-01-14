process SOURMASH_SIG_TO_ZIP {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::sourmash=4.9.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.9.4--hdfd78af_0':
        'biocontainers/sourmash:4.9.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(sig_file)

    output:
    tuple val(meta), path("*.sig.zip"), emit: sig_zip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash sig cat \\
        ${sig_file} \\
        -o ${prefix}.sig.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sig.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
