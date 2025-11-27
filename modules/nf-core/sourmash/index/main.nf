process SOURMASH_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"

    input:
    tuple val(meta), path(signatures)
    val(ksize)

    output:
    tuple val(meta), path("*.sbt.zip"), emit: signature_index
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash index \\
        --ksize ${ksize} \\
        $args \\
        '${prefix}.sbt.zip' \\
        $signatures

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.sbt.zip"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
