process SOURMASH_FASTMULTIGATHER {
    tag "${manysketch_zip.baseName}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"

    input:
    tuple val(meta), path(manysketch_zip)
    path sourmash_database

    output:
    tuple val(meta), path("fastmultigather/*_sourmash_gather.csv"), emit: gather_csv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${manysketch_zip.baseName}"
    def args = task.ext.args ?: ''
    def threshold = params.sourmash_thresholdbp ? "--threshold-bp=${params.sourmash_thresholdbp}" : '--threshold-bp=50000'
    """
    sourmash scripts fastmultigather \\
        ${manysketch_zip} \\
        ${sourmash_database} \\
        -c ${task.cpus} \\
        -k ${params.sourmash_ksize} \\
        -o ${prefix}_sourmash_gather_withrocksdb.csv \\
        ${threshold} \\
        $args 2>&1 | tee fastmultigather.log

    mkdir -p fastmultigather
    mv ${prefix}_sourmash_gather_withrocksdb.csv fastmultigather/${prefix}_sourmash_gather.csv

    # Extract no match bins
    if [ -f fastmultigather.log ]; then
        grep "No matches to" fastmultigather.log | awk -F"'" '{print \$2}' > fastmultigather/no_match_bins.txt || touch fastmultigather/no_match_bins.txt
        mv fastmultigather.log fastmultigather/
    fi

    # Move any prefetch files if they exist
    find . -name "*.prefetch.csv" -exec mv {} fastmultigather/ \\; || true

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
    def prefix = task.ext.prefix ?: "fastmultigather_results"
    """
    mkdir -p fastmultigather
    touch fastmultigather/${prefix}_sourmash_gather.csv
    touch fastmultigather/no_match_bins.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: "stub_version"
        sed: "stub_version"
    END_VERSIONS
    """
}