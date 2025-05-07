process PROCESS_RESULTS {
    tag "Processing results for ${meta.id}"
    label 'process_medium'

    conda "${projectDir}/env/read_based.yaml"

    publishDir "${params.outdir}/processed_results",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(gather_csv)
    tuple val(meta), path(yacht_xlsx)

    output:
    tuple val(meta), path("final_results"), emit: final_results
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: 'results'
    """
    mkdir -p final_results
    mkdir -p temp_text_files
    mkdir -p temp_yacht_results

    cp ${gather_csv} temp_text_files/sourmash_gather_withrocksdb.csv
    cp ${yacht_xlsx} temp_yacht_results/

    python ${projectDir}/bin/py_scripts/1_process_sourmash.py \\
        -i temp_text_files/sourmash_gather_withrocksdb.csv \\
        -o temp_text_files/sourmash_mergedresults.csv

    python ${projectDir}/bin/py_scripts/2_process_yacht.py \\
        -d temp_yacht_results

    python ${projectDir}/bin/py_scripts/3_concat_yacht.py \\
        -d temp_yacht_results

    python ${projectDir}/bin/py_scripts/4_merged_yacht_sourmash.py \\
        -d1 final_results \\
        -d2 temp_text_files \\
        -d3 temp_yacht_results

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