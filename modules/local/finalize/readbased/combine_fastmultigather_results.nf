process COMBINE_FASTMULTIGATHER_RESULTS {
    tag "${meta.id}"
    label 'process_low'

    conda "coreutils=9.1"

    input:
    tuple val(meta), path(csv_files)

    output:
    tuple val(meta), path("*_combined*.csv"), emit: combined_csv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = meta.file_type ?: 'gather'
    
    // Determine output filename based on file type
    def output_name = file_type == 'taxannotate' ? "${prefix}_combined_sourmash_gather.with-lineages.csv" :
                      file_type == 'taxgenome' ? "${prefix}_combined_genome_classification.csv" :
                      "${prefix}_combined_sourmash_gather.csv"
    
    """
    #!/bin/bash
    
    # Count number of input files
    num_files=\$(echo ${csv_files} | wc -w)
    
    echo "Combining \${num_files} file(s) for sample: ${meta.id} (type: ${file_type})"
    
    # Decompress .gz files if present
    for file in ${csv_files}; do
        if [[ "\${file}" == *.gz ]]; then
            echo "Decompressing: \${file}"
            gunzip -f "\${file}"
        fi
    done
    
    # Get all CSV files (now decompressed)
    csv_list=\$(ls *.csv 2>/dev/null || echo "")
    
    if [ -z "\${csv_list}" ]; then
        echo "ERROR: No CSV files found after decompression"
        exit 1
    fi
    
    # Get the first CSV file to extract header
    first_csv=\$(echo \${csv_list} | cut -d' ' -f1)
    
    # Write header from first file
    head -n 1 "\${first_csv}" > ${output_name}
    
    # Append all data rows (skip header) from all CSV files
    for csv_file in \${csv_list}; do
        echo "Processing: \${csv_file}"
        tail -n +2 "\${csv_file}" >> ${output_name}
    done
    
    # Count total rows
    total_rows=\$(tail -n +2 ${output_name} | wc -l)
    echo "Total rows combined: \${total_rows}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n 1 | cut -d ' ' -f 4)
        gunzip: \$(gunzip --version 2>&1 | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = meta.file_type ?: 'gather'
    def output_name = file_type == 'taxannotate' ? "${prefix}_combined_sourmash_gather.with-lineages.csv" :
                      file_type == 'taxgenome' ? "${prefix}_combined_genome_classification.csv" :
                      "${prefix}_combined_sourmash_gather.csv"
    """
    echo "query_filename,query_name,query_md5" > ${output_name}
    echo "test_bin,test,abc123" >> ${output_name}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub_version"
        gunzip: "stub_version"
    END_VERSIONS
    """
}