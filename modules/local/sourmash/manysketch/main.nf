process SOURMASH_MANYSKETCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.sketches_dir}", mode: params.publish_dir_mode, enabled: params.enable_copysketch

    input:
    tuple val(meta), path(manysketch_csv)

    output:
    tuple val(meta), path("*.manysketch.zip"), emit: sketch_zip_file
    path "versions.yml"                      , emit: versions

    script:
    def prefix = task.ext.prefix ?: "manysketch_results"
    def args_sourmash = task.ext.args_sourmash ?: ''
    """
    # Run sourmash manysketch using the provided CSV file
    sourmash scripts manysketch ${manysketch_csv} \\
        -o ${prefix}.batch.manysketch.zip \\
        -c ${task.cpus} \\
        -p ${params.sourmash_profile} \\
        $args_sourmash

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //g; s/ version//g' | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch batch.manysketch.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: "stub_version"
    END_VERSIONS
    """
}