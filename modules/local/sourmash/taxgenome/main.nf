process SOURMASH_TAXGENOME {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"

    input:
    tuple val(meta), path(taxannotate_csv_gz)
    path sourmash_taxonomy_csv

    output:
    tuple val(meta), path("taxgenome/*.csv"), emit: genome_classification
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: "taxgenome_results"
    def args = task.ext.args ?: ''
    def output_format = params.sourmash_output_format ?: 'csv_summary'
    """
    # Decompress the input file first
    gunzip -c ${taxannotate_csv_gz} > gather_with_lineages.csv
    
    mkdir -p taxgenome
    sourmash tax genome \\
        --gather-csv gather_with_lineages.csv \\
        --taxonomy ${sourmash_taxonomy_csv} \\
        --output-dir taxgenome \\
        --output-format ${output_format} \\
        --output-base ${prefix} \\
        $args
    
    # Rename output file to match expected pattern
    if [ -f "taxgenome/${prefix}.${output_format}.csv" ]; then
        mv "taxgenome/${prefix}.${output_format}.csv" "taxgenome/${prefix}_genome_classification.csv"
    elif [ -f "taxgenome/${prefix}.csv" ]; then
        mv "taxgenome/${prefix}.csv" "taxgenome/${prefix}_genome_classification.csv"
    fi

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