process SOURMASH_MANYSKETCH {
    tag "Running sourmash manysketch on batch"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.sketches_dir}", mode: params.publish_dir_mode, enabled: params.enable_copysketch

    input:
    path manysketch_csv

    output:
    path "manysketch_output"                         , emit: manysketch_dir
    path "manysketch_output/zip_files"              , emit: zip_files_dir
    path "batch.manysketch.zip"                     , emit: sketch_zip_file
    path "versions.yml"                             , emit: versions

    script:
    def args_sourmash = task.ext.args_sourmash ?: '' // Args for sourmash manysketch
    def args_python = task.ext.args_python ?: ''     // Args for the python script
    """
    # Run sourmash manysketch using the provided CSV file
    sourmash scripts manysketch ${manysketch_csv} \\
        -o batch.manysketch.zip \\
        -c ${task.cpus} \\
        -p dna,k=31,k=51,scaled=1000,abund \\
        $args_sourmash

    # Create directories
    mkdir -p manysketch_output/zip_files

    # Unzip into manysketch_output
    unzip -o batch.manysketch.zip -d manysketch_output/

    # Run the Python script on the unzipped output
    python ${projectDir}/bin/py_scripts/read_based/sourmash_extract_zip_with_ksize.py -d manysketch_output $args_python

    # Move any generated zip files to zip_files directory
    mv manysketch_output/*_${params.ksize}.sig.zip manysketch_output/zip_files/ 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //g; s/ version//g' | head -n 1)
        python: \$(python --version 2>&1 | sed 's/Python //g')
        unzip: \$(unzip -v 2>&1 | head -n 1 | sed 's/UnZip //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p manysketch_output/zip_files
    touch batch.manysketch.zip
    touch manysketch_output/zip_files/test_${params.ksize}.sig.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: "stub_version"
        python: "stub_version"
        unzip: "stub_version"
    END_VERSIONS
    """
}