process SOURMASH_TAXMETAGENOME {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"
    
    input:
    tuple val(meta), path(taxannotate_csv_gz)
    path sourmash_taxonomy_csv

    output:
    tuple val(meta), path("taxmetagenome/*.csv"), emit: genome_classification
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: "taxmetagenome_results"
    def args = task.ext.args ?: ''
    def output_format = params.sourmash_output_format ?: 'csv_summary'
    """
    # Decompress the input file first
    gunzip -c ${taxannotate_csv_gz} > gather_with_lineages.csv
    
    mkdir -p taxmetagenome
    sourmash tax metagenome \\
        --gather-csv gather_with_lineages.csv \\
        --taxonomy ${sourmash_taxonomy_csv} \\
        --output-dir taxmetagenome \\
        --output-format ${output_format} \\
        --output-base ${prefix} \\
        $args
    
    # Rename output file to match expected pattern
    if [ -f "taxmetagenome/${prefix}.${output_format}.csv" ]; then
        mv "taxmetagenome/${prefix}.${output_format}.csv" "taxmetagenome/${prefix}_metagenome_classification.csv"
    elif [ -f "taxmetagenome/${prefix}.csv" ]; then
        mv "taxmetagenome/${prefix}.csv" "taxmetagenome/${prefix}_metagenome_classification.csv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //g; s/ version//g' | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}"
    """
    mkdir -p taxmetagenome
    touch taxmetagenome/${prefix}_metagenome_classification.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: "stub_version"
    END_VERSIONS
    """
}