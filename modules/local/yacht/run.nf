process YACHT_RUN {
    tag "Running YACHT on batch samples"
    label 'process_high'

    conda "${projectDir}/env/read_based.yaml"

    publishDir "${params.yacht_results_dir}",
        mode: params.publish_dir_mode,
        enabled: params.enable_copyintermediate,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(zip_files_dir)  // Directory containing zip files
    path yacht_database_json

    output:
    tuple val(meta), path("yacht_results/*.xlsx"), emit: yacht_xlsx
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: 'yacht'
    def args = task.ext.args ?: '' // Added for custom arguments to yacht run
    """
    mkdir -p yacht_results
    
    # Copy signature files to current directory
    cp ${zip_files_dir}/*_${params.ksize}.sig.zip .
    
    # Process each signature file
    for sig_zip in *_${params.ksize}.sig.zip; do
        sample_name=\$(basename \$sig_zip _${params.ksize}.sig.zip)
        
        yacht run \\
            --json ${yacht_database_json} \\
            --sample_file \$sig_zip \\
            --num_threads ${task.cpus} \\
            --significance 0.99 \\
            --min_coverage_list 1 0.5 0.1 0.05 0.01 \\
            --out yacht_results/\${sample_name}_yacht.xlsx \\
            $args
    done
    
    # Cleanup
    rm -f *_${params.ksize}.sig.zip
    find . -type d -name "*intermediate_files" -exec rm -rf {} +

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yacht: \$(yacht --version 2>&1 | head -n 1)
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