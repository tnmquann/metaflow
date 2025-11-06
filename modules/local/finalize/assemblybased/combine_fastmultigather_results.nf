process COMBINE_FASTMULTIGATHER_RESULTS {
    tag "Combining ${gather_csvs.size()} result files for ${meta.id}"
    label 'process_low'

    conda "coreutils=9.1"

    input:
    tuple val(meta), path(gather_csvs)

    output:
    tuple val(meta), path("*_combined_sourmash_gather.csv"), emit: gather_csv
    path "versions.yml", emit: versions

    script:
    // Use meta.id as prefix since we're combining across all bins for this sample
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    
    # Count number of input CSV files
    num_files=\$(echo ${gather_csvs} | wc -w)
    
    echo "Combining \${num_files} CSV file(s) for sample: ${meta.id}"
    
    # Get the first CSV file to extract header
    first_csv=\$(echo ${gather_csvs} | cut -d' ' -f1)
    
    # Write header from first file
    head -n 1 "\${first_csv}" > ${prefix}_combined_sourmash_gather.csv
    
    # Append all data rows (skip header) from all CSV files
    for csv_file in ${gather_csvs}; do
        echo "Processing: \${csv_file}"
        tail -n +2 "\${csv_file}" >> ${prefix}_combined_sourmash_gather.csv
    done
    
    # Count total matches
    total_matches=\$(tail -n +2 ${prefix}_combined_sourmash_gather.csv | wc -l)
    echo "Total matches combined: \${total_matches}"
    
    # Optional: Sort by f_unique_weighted (column 16) in descending order
    if [ \${total_matches} -gt 0 ]; then
        head -n 1 ${prefix}_combined_sourmash_gather.csv > temp_header.csv
        tail -n +2 ${prefix}_combined_sourmash_gather.csv | sort -t',' -k16 -gr > temp_sorted.csv 2>/dev/null || tail -n +2 ${prefix}_combined_sourmash_gather.csv > temp_sorted.csv
        cat temp_header.csv temp_sorted.csv > ${prefix}_combined_sourmash_gather.csv
        rm -f temp_header.csv temp_sorted.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n 1 | cut -d ' ' -f 4)
        sort: \$(sort --version 2>&1 | head -n 1 | cut -d ' ' -f 4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub combined CSV
    echo "query_filename,query_name,query_md5,query_bp,ksize,moltype,scaled,query_n_hashes,filename,name,md5,f_match_orig,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund" > ${prefix}_combined_sourmash_gather.csv
    echo "test_bin,test,abc123,50000,51,DNA,1000,50,db.zip,match,def456,0.9,0.1,0.9,0.2,0.8,10.0,8.0,3.0" >> ${prefix}_combined_sourmash_gather.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub_version"
        sort: "stub_version"
    END_VERSIONS
    """
}