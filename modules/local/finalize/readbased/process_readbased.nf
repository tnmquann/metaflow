process PROCESS_READBASED_RESULTS {
    tag "Processing results for ${meta.id}"
    label 'process_medium'

    conda "conda-forge::pandas=2.2.3 openpyxl=3.1.5"

    publishDir "${params.processed_results_dir}",
        mode: params.publish_dir_mode,
        pattern: "final_results"

    publishDir "${params.yacht_results_dir}",
        mode: params.publish_dir_mode,
        pattern: "temp_yacht_results",
        saveAs: { 'yacht_results' }

    input:
    tuple val(meta), path(gather_csv)
    tuple val(meta), path(yacht_xlsx)

    output:
    tuple val(meta), path("final_results"), emit: final_results
    path "versions.yml", emit: versions
    path "temp_yacht_results", optional: true

    script:
    def prefix = task.ext.prefix ?: meta.id ?: 'results'
    def py_args_script1 = task.ext.py_args_script1 ?: ''
    def py_args_script2 = task.ext.py_args_script2 ?: ''
    def py_args_script3 = task.ext.py_args_script3 ?: ''
    def py_args_script4 = task.ext.py_args_script4 ?: ''
    """
    mkdir -p final_results
    mkdir -p temp_text_files
    mkdir -p temp_yacht_results

    cp ${gather_csv} temp_text_files/sourmash_gather_withrocksdb.csv
    cp ${yacht_xlsx} temp_yacht_results/

    python ${projectDir}/bin/py_scripts/read_based/process_sourmash_results.py \\
        -i temp_text_files/sourmash_gather_withrocksdb.csv \\
        -o temp_text_files/sourmash_mergedresults.csv \\
        $py_args_script1

    python ${projectDir}/bin/py_scripts/read_based/process_yacht_results.py \\
        -d temp_yacht_results \\
        $py_args_script2

    python ${projectDir}/bin/py_scripts/read_based/concat_yacht_xlsx.py \\
        -d temp_yacht_results \\
        $py_args_script3

    python ${projectDir}/bin/py_scripts/read_based/merged_yacht_sourmash.py \\
        -d1 final_results \\
        -d2 temp_text_files \\
        -d3 temp_yacht_results \\
        $py_args_script4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id ?: 'results'
    """
    mkdir -p final_results
    touch final_results/summary_${prefix}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub_version"
    END_VERSIONS
    """
}