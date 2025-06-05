process SOURMASH_FASTMULTIGATHER {
    tag "Running sourmash fastmultigather on ${manysketch_zip}"
    label 'process_high'

    conda "${projectDir}/env/read_based.yaml"

    publishDir "${params.fastmultigather_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(manysketch_zip)
    path sourmash_database

    output:
    tuple val(meta), path("fastmultigather/*_sourmash_gather.csv"), emit: gather_csv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id ?: "fastmultigather_results"
    """
    sourmash scripts fastmultigather \\
        ${manysketch_zip} \\
        ${sourmash_database} \\
        -c ${task.cpus} \\
        -k ${params.ksize} \\
        -o ${prefix}_sourmash_gather_withrocksdb.csv

    mkdir -p fastmultigather
    mv ${prefix}_sourmash_gather_withrocksdb.csv fastmultigather/${prefix}_sourmash_gather.csv

    cd fastmultigather
    sed -i 's/match_filename/filename/g' ${prefix}_sourmash_gather.csv
    sed -i 's/match_name/name/g' ${prefix}_sourmash_gather.csv
    sed -i 's/match_md5/md5/g' ${prefix}_sourmash_gather.csv
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //g; s/ version//g' | head -n 1)
        sed: \$(sed --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id ?: "fastmultigather_results"
    """
    mkdir -p fastmultigather
    touch fastmultigather/${prefix}_sourmash_gather.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: "stub_version"
        sed: "stub_version"
    END_VERSIONS
    """
}