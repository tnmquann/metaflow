process PROCESS_RESULTS {
    tag "Processing final results"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path text_files
    path yacht_results 
    path merged_seq
    path manysketch
    path sketch_zip

    output:
    path "final_results"

    script:
    """
    # Create final_results directory
    mkdir -p final_results

    # Process the data first
    python ${params.workDir}/bin/py_scripts/1_process_sourmash.py -i ${text_files}/sourmash_gather_withrocksdb.csv -o ${text_files}/sourmash_mergedresults.csv
    python ${params.workDir}/bin/py_scripts/2_process_yacht.py -d ${yacht_results}
    python ${params.workDir}/bin/py_scripts/3_concat_yacht.py -d ${yacht_results}
    python ${params.workDir}/bin/py_scripts/4_merged_yacht_sourmash.py -d1 final_results -d2 ${text_files} -d3 ${yacht_results}

    # Create all target directories in base_dir
    mkdir -p "${params.base_dir}/final_results"
    if [ "${params.enable_copymergedseqs}" = "true" ]; then
        mkdir -p "${params.base_dir}/merged_seq"
    fi
    if [ "${params.enable_copysketch}" = "true" ]; then
        mkdir -p "${params.base_dir}/sketches"
    fi
    if [ "${params.enable_copyintermediate}" = "true" ]; then
        mkdir -p "${params.base_dir}/intermediate_files"
    fi

    # Copy files using cp -rL to resolve symlinks
    cp -rL final_results/* "${params.base_dir}/final_results/"

    # Copy other directories if enabled, using cp -rL to resolve symlinks
    if [ "${params.enable_copymergedseqs}" = "true" ]; then
        cp -rL "${merged_seq}" "${params.base_dir}/merged_seq/"
    fi

    if [ "${params.enable_copysketch}" = "true" ]; then
        cp -rL "${manysketch}" "${params.base_dir}/sketches/"
        cp -L "${sketch_zip}" "${params.base_dir}/sketches/"
    fi

    if [ "${params.enable_copyintermediate}" = "true" ]; then
        cp -rL "${text_files}" "${params.base_dir}/intermediate_files/"
        cp -rL "${yacht_results}" "${params.base_dir}/intermediate_files/"
    fi
    """
}