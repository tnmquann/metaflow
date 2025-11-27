process SOURMASH_TAXANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"

    input:
    tuple val(meta), path(gather_results)
    path(taxonomy)

    output:
    tuple val(meta), path("*.with-lineages.csv.gz"), emit: result
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash \\
        tax annotate \\
        $args \\
        --gather-csv ${gather_results} \\
        --taxonomy ${taxonomy} \\
        --output-dir "."

    ## Compress output
    gzip --no-name *.with-lineages.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > test.with-lineages.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
