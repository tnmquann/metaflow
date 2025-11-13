process EXTRACT_SOURMASH_SINGLESKETCHES {
    tag "${meta.id}"
    label 'process_low'

    conda "pandas=2.2.3"

    input:
    tuple val(meta), path(sketch_zip)
    val ksize

    output:
    tuple val(meta), path("manysketch_output")                  , emit: manysketch_dir, optional: true
    tuple val(meta), path("manysketch_output/zip_files/*.sig.zip"), emit: sig_zip_files, optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args_python = task.ext.args_python ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create directories
    mkdir -p manysketch_output/zip_files

    # Unzip into manysketch_output
    unzip -o ${sketch_zip} -d manysketch_output/

    # Run the Python script on the unzipped output
    python ${projectDir}/bin/py_scripts/read_based/sourmash_extract_zip_with_ksize.py \\
        -d manysketch_output \\
        $args_python

    # Move any generated zip files to zip_files directory
    mv manysketch_output/*_${ksize}.sig.zip manysketch_output/zip_files/ 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        unzip: \$(unzip -v 2>&1 | head -n 1 | sed 's/UnZip //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p manysketch_output/zip_files
    touch manysketch_output/zip_files/${prefix}_test_${ksize}.sig.zip
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub_version"
        pandas: "stub_version"
        unzip: "stub_version"
    END_VERSIONS
    """
}
