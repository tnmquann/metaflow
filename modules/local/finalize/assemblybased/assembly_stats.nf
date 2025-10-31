process ASSEMBLY_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.3 conda-forge::biopython=1.86"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")

    output:
    tuple val(meta), path("*.assembly_stats.csv"), emit: asm_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/py_scripts/assembly_based/assembly_stats.py \\
         --fasta_dir input_bins/ \\
         --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python ${projectDir}/bin/py_scripts/assembly_based/assembly_stats.py: \$(assembly_stats.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}
