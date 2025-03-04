process YACHT {
    tag "Running YACHT"
    publishDir "${params.base_dir}/intermediate_files",
        mode: 'copy',
        enabled: params.enable_copyintermediate
    
    input:
    path zip_files
    val yacht_database  // changed from 'path yacht_database'

    output:
    path "yacht_results"

    script:
    """
    mkdir -p yacht_results
    cp ${zip_files}/*_${params.ksize}.sig.zip .
    
    for i in *_${params.ksize}.sig.zip; do
        k=\$(basename \$i _${params.ksize}.sig.zip)
        yacht run --json ${yacht_database} --sample_file \$i --num_threads ${task.cpus} --significance 0.99 --min_coverage_list 0.5 0.1 0.05 --out yacht_results/\${k}_yacht.xlsx
    done
    
    rm *_${params.ksize}.sig.zip
    find . -type d -name "*intermediate_files" -exec rm -rf {} +
    """
}