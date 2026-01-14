process YACHT_RUN_SINGLESKETCH {
    tag "Running YACHT on ${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/sourmash_yacht_pandas_sourmash_plugin_branchwater:f5cacb08e855b10c"

    publishDir "${params.yacht_results_dir}",
        mode: params.publish_dir_mode,
        enabled: params.enable_copyintermediate,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(sig_file)  // Single signature file
    path yacht_database_json

    output:
    tuple val(meta), path("yacht_results/${meta.id}_yacht.xlsx"), emit: yacht_xlsx
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: 'yacht'
    def args = task.ext.args ?: ''
    """
    mkdir -p yacht_results
    
    # Run YACHT directly on the signature file
    yacht run \\
        --json ${yacht_database_json} \\
        --sample_file ${sig_file} \\
        --num_threads ${task.cpus} \\
        --significance 0.99 \\
        --min_coverage_list 1 0.5 0.1 0.05 0.01 \\
        --out yacht_results/${prefix}_yacht.xlsx \\
        $args
    
    # Cleanup intermediate files
    find . -type d -name "*intermediate_files" -exec rm -rf {} +

    # Correct version extraction for YACHT
    YACHT_VERSION=\$(yacht --version | grep -oP '(?<=version )[^ ]+(?= )')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yacht: "\$YACHT_VERSION"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id ?: 'yacht'
    """
    mkdir -p yacht_results
    touch yacht_results/${prefix}_yacht.xlsx
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yacht: "stub_version"
    END_VERSIONS
    """
}
