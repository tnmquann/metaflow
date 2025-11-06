process SOURMASH_TAXGENOME {
    tag "Running sourmash tax genome on ${gather_csv.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gather_csv)
    path sourmash_taxonomy_csv

    output:
    tuple val(meta), path("taxgenome/*.csv"), emit: genome_classification
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: "taxgenome_results"
    def args = task.ext.args ?: ''
    def output_format = params.sourmash_output_format ?: 'csv_summary'
    """
    mkdir -p taxgenome
    sourmash tax genome \\
        --gather-csv ${gather_csv} \\
        --taxonomy ${sourmash_taxonomy_csv} \\
        --output-dir taxgenome \\
        --output-format ${output_format} \\
        --output-base ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //g; s/ version//g' | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}"
    """
    mkdir -p taxgenome
    touch taxgenome/${prefix}_genome_classification.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: "stub_version"
    END_VERSIONS
    """
}