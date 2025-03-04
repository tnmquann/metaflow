process SOURMASH_MANYSKETCH {
    tag "Running sourmash manysketch"
    publishDir "${params.base_dir}/sketches",
        mode: 'copy',
        enabled: params.enable_copysketch

    input:
    path merged_seq

    output:
    path "manysketch", emit: manysketch_dir
    path "manysketch/zip_files", emit: zip_files
    path "manysketch.zip", emit: sketch_zip

    script:
    """
    echo "name,genome_filename,protein_filename" > manysketch.csv
    for i in ${merged_seq}/*.fastq.gz; do
        k=\$(basename \$i .fastq.gz)
        echo "\$k,\$i," >> manysketch.csv
    done

    sourmash scripts manysketch manysketch.csv -o manysketch.zip -c ${task.cpus} -p dna,k=31,k=51,scaled=1000,abund

    unzip -o manysketch.zip -d manysketch
    python ${params.workDir}/bin/py_scripts/0_zip_with_ksize_ok.py -d manysketch
    """
}