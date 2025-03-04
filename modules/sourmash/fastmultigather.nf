process SOURMASH_FASTMULTIGATHER {
    tag "Running sourmash fastmultigather"
    publishDir "${params.base_dir}/intermediate_files", 
        mode: 'copy',
        enabled: params.enable_copyintermediate
    
    input:
    path manysketch_zip
    path sourmash_database
    path merged_seq

    output:
    path "text_files"

    script:
    """
    sourmash scripts fastmultigather ${manysketch_zip} ${sourmash_database} -c ${task.cpus} -k 51 -o sourmash_gather_withrocksdb.csv
    
    mkdir -p text_files
    mv sourmash_gather_withrocksdb.csv text_files/
    
    cd text_files
    sed -i 's/match_filename/filename/g' sourmash_gather_withrocksdb.csv
    sed -i 's/match_name/name/g' sourmash_gather_withrocksdb.csv
    sed -i 's/match_md5/md5/g' sourmash_gather_withrocksdb.csv
    """
}