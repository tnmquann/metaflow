process RGI_PREPARECARDDB {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1':
        'biocontainers/rgi:6.0.3--pyha8f3691_1' }"

    output:
    path("rgi_db") , emit: db
    env RGI_VERSION                 , emit: tool_version
    env DB_VERSION                  , emit: db_version
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Download and extract CARD data
    mkdir -p ./rgi_db
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data -C ./rgi_db/
    cd ./rgi_db

    # Run card annotation
    rgi card_annotation \\
        -i ./card.json > card_annotation.log 2>&1

    # Extract database version
    DB_VERSION=\$(ls card_database_*_all.fasta | sed "s/card_database_v\\([0-9].*[0-9]\\).*/\\1/")

    # Download and extract wildcard data
    wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
    mkdir -p wildcard
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    gunzip wildcard/*.gz

    # Run wildcard annotation
    rgi wildcard_annotation \\
        -i wildcard \\
        --card_json ./card.json \\
        -v \${DB_VERSION} > wildcard_annotation.log 2>&1

    # Load the database
    rgi load \\
        --card_json ./card.json \\
        --wildcard_annotation wildcard_database_v\${DB_VERSION}.fasta \\
        --wildcard_annotation_all_models wildcard_database_v\${DB_VERSION}_all.fasta \\
        --wildcard_index ./wildcard/index-for-model-sequences.txt \\
        --wildcard_version \${DB_VERSION} \\
        --card_annotation card_database_v\${DB_VERSION}.fasta \\
        --card_annotation_all_models card_database_v\${DB_VERSION}_all.fasta \\
        --debug \\
        --local \\
        ${args}

    cd ..
    RGI_VERSION=\$(rgi main --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ./rgi_db
    cd ./rgi_db

    # Create minimal CARD files expected by the real process
    touch card.json
    touch card_database_v4.0.1.fasta
    touch card_database_v4.0.1_all.fasta

    # Create wildcard dir and files used by the real process
    mkdir -p wildcard
    touch wildcard/wildcard_database_v4.0.1.fasta
    touch wildcard/wildcard_database_v4.0.1_all.fasta
    mkdir -p wildcard
    touch wildcard/index-for-model-sequences.txt

    # Do not move/copy files out of rgi_db — rgi_db is the emitted output
    cd ..

    # Try to get real rgi version if available, otherwise fallback
    RGI_VERSION=\$(rgi main --version 2>/dev/null || echo "rgi-stub")
    DB_VERSION=stub_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """
}