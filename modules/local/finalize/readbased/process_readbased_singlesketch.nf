process PROCESS_READBASED_RESULTS_SINGLESKETCH {
    tag "Processing results for ${meta.id}"
    label 'process_medium'

    conda "conda-forge::pandas=2.2.3 openpyxl=3.1.5"

    input:
    tuple val(meta), path(gather_csv)
    tuple val(meta2), path(yacht_xlsx)

    output:
    tuple val(meta), path("*_processedresult"), emit: final_results
    path "versions.yml", emit: versions
    path "yacht_processed", optional: true

    script:
    def prefix = task.ext.prefix ?: meta.id ?: 'sample'
    def min_coverage = task.ext.min_coverage ?: 0.05
    def py_args = task.ext.py_args ?: ''
    """
    mkdir -p ${prefix}
    mkdir -p yacht_processed

    # Decompress the gather CSV file if needed
    if [[ ${gather_csv} == *.gz ]]; then
        gunzip -c ${gather_csv} > gather_input.csv
        GATHER_INPUT="gather_input.csv"
    else
        GATHER_INPUT="${gather_csv}"
    fi

    # Run the single-sample processing script
    python ${projectDir}/bin/py_scripts/read_based/process_singlesketch_results.py \\
        -s ${prefix} \\
        -g \$GATHER_INPUT \\
        -y ${yacht_xlsx} \\
        -o ${prefix}_processedresult \\
        --remove_unclassified \\
        --min_coverage ${min_coverage} \\
        $py_args

    # Copy YACHT xlsx to yacht_processed for publishing
    cp ${yacht_xlsx} yacht_processed/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        openpyxl: \$(python -c "import openpyxl; print(openpyxl.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id ?: 'sample'
    """
    mkdir -p ${prefix}
    mkdir -p yacht_processed
    touch ${prefix}/${prefix}_merged_sourmash_yacht.csv
    touch yacht_processed/${prefix}_yacht.xlsx
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub_version"
    END_VERSIONS
    """
}