process RGI_PREPARECARDDB {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/rgi_aria2:a94e11e18f7348f4':
        'biocontainers/rgi:6.0.5--pyh05cac1d_0' }"

    input:
    path card, stageAs: 'input_card/*'  // Optional input

    output:
    path("rgi_db"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def download_threads = 8

    """
    # Create output directory
    mkdir -p ./rgi_db
    DB_VERSION=""
    RGI_VERSION=""

    download_file() {
        local output_file="\$1"
        local url="\$2"

        if command -v aria2c >/dev/null 2>&1; then
            aria2c \\
                --allow-overwrite=true \\
                --continue=true \\
                --file-allocation=none \\
                --max-connection-per-server=${download_threads} \\
                --split=${download_threads} \\
                --min-split-size=1M \\
                --dir=. \\
                --out="\${output_file}" \\
                "\${url}"
        else
            wget -O "\${output_file}" "\${url}"
        fi
    }

    # Check if input was provided
    if [ -d "input_card" ] && [ "\$(ls -A input_card 2>/dev/null)" ]; then
        echo "Using provided CARD data..."
        
        # Copy input data to working directory
        cp -r input_card/* ./rgi_db/
        cd ./rgi_db
        
        # Verify required files exist
        if [ ! -f "card.json" ]; then
            echo "ERROR: card.json not found in input data" >&2
            exit 1
        fi
        
        # Extract database version from existing files
        DB_VERSION=\$(find . -name "card_database_*_all.fasta" -exec basename {} \\; | \\
                    sed 's/card_database_v\\([0-9.]*\\).*/\\1/' | head -1)
        if [ -z "\$DB_VERSION" ]; then
            echo "ERROR: Could not determine database version from input files" >&2
            exit 1
        fi
        
        # Verify wildcard files exist
        if [ ! -d "wildcard" ] || [ ! -f "wildcard/index-for-model-sequences.txt" ]; then
            echo "ERROR: Required wildcard files not found in input data" >&2
            exit 1
        fi
        
        # Only run rgi load with existing data
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
            ${args} || {
            echo "ERROR: RGI database loading failed" >&2
            exit 1
        }
        
    else
        echo "No input provided, downloading and processing CARD data..."
        
        # Download and extract CARD data
        download_file card_data.tar.bz2 https://card.mcmaster.ca/latest/data || {
            echo "ERROR: Failed to download CARD data" >&2
            exit 1
        }
        tar -xjf card_data.tar.bz2 -C ./rgi_db/ || {
            echo "ERROR: Failed to extract CARD data" >&2
            exit 1
        }
        rm -f card_data.tar.bz2
        cd ./rgi_db

        # Verify card.json exists
        if [ ! -f "card.json" ]; then
            echo "ERROR: card.json not found in downloaded data" >&2
            exit 1
        fi

        # Run card annotation
        rgi card_annotation \\
            -i ./card.json > card_annotation.log 2>&1 || {
            echo "ERROR: card_annotation failed" >&2
            exit 1
        }

        # Extract database version
        DB_VERSION=\$(find . -name "card_database_*_all.fasta" -exec basename {} \\; | \\
                    sed 's/card_database_v\\([0-9.]*\\).*/\\1/' | head -1)
        if [ -z "\$DB_VERSION" ]; then
            echo "ERROR: Could not determine database version" >&2
            exit 1
        fi

        # Download and extract wildcard data
        download_file wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants || {
            echo "ERROR: Failed to download wildcard data" >&2
            exit 1
        }
        mkdir -p wildcard
        tar -xjf wildcard_data.tar.bz2 -C wildcard || {
            echo "ERROR: Failed to extract wildcard data" >&2
            exit 1
        }
        rm -f wildcard_data.tar.bz2

        # Extract wildcard .gz files if they exist
        if [ -d "wildcard" ] && [ "\$(ls -A wildcard/*.gz 2>/dev/null)" ]; then
            gunzip wildcard/*.gz
        fi

        # Run wildcard annotation
        rgi wildcard_annotation \\
            -i wildcard \\
            --card_json ./card.json \\
            -v \${DB_VERSION} > wildcard_annotation.log 2>&1 || {
            echo "ERROR: wildcard_annotation failed" >&2
            exit 1
        }

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
            ${args} || {
            echo "ERROR: RGI database loading failed" >&2
            exit 1
        }
    fi

    # Get tool version
    cd ..
    RGI_VERSION=\$(rgi main --version 2>/dev/null || echo "unknown")
    export DB_VERSION
    export RGI_VERSION

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \${RGI_VERSION}
        rgi-database: \${DB_VERSION}
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

    # Create wildcard directory and files
    mkdir -p wildcard
    touch wildcard/wildcard_database_v4.0.1.fasta
    touch wildcard/wildcard_database_v4.0.1_all.fasta
    touch wildcard/index-for-model-sequences.txt

    # Create localDB directory (simulating what rgi load --local would create)
    mkdir -p localDB
    touch localDB/card.json

    cd ..

    # Get version with fallback
    RGI_VERSION=\$(rgi main --version 2>/dev/null || echo "rgi-stub")
    DB_VERSION="4.0.1"
    export DB_VERSION
    export RGI_VERSION

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \${RGI_VERSION}
        rgi-database: \${DB_VERSION}
    END_VERSIONS
    """
}
