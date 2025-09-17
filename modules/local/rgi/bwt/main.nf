process RGI_BWT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1':
        'biocontainers/rgi:6.0.3--pyha8f3691_1' }"

    input:
    tuple val(meta), path(reads1), path(reads2)  // meta, forward_reads, reverse_reads
    path db                                       // RGI database directory

    output:
    tuple val(meta), path("${meta.id}"), emit: outdir        // Changed from ${prefix}/ to ${meta.id}
    tuple val(meta), path("${meta.id}/*.json"), emit: json, optional: true
    tuple val(meta), path("${meta.id}/*.txt"), emit: tsv, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Determine DB version from prepared DB folder (look for _all.fasta files)
    DB_VERSION=\$(find ${db}/ -name "card_database_*_all.fasta" -exec basename {} \\; | \\
                  sed 's/card_database_v\\([0-9.]*\\)_all.*/\\1/' | head -1)
    if [ -z "\$DB_VERSION" ]; then
        echo "WARNING: Could not determine database version from _all.fasta files, trying regular fasta files"
        DB_VERSION=\$(find ${db}/ -name "card_database_*.fasta" -exec basename {} \\; | \\
                      sed 's/card_database_v\\([0-9.]*\\).*/\\1/' | head -1 || echo "unknown")
    fi

    # Check if localDB exists before creating symlink
    if [ -d "${db}/localDB" ]; then
        # Create symlink with standard name that RGI expects
        ln -s \$(realpath ${db}/localDB) ./localDB
    else
        echo "ERROR: localDB not found in ${db}" >&2
        exit 1
    fi

    # Create output directory
    mkdir -p ${prefix}
    
    # Set default RGI bwt options (can be overridden by args)
    DEFAULT_ARGS="--local -a kma --clean"
    
    # Run RGI BWT
    rgi bwt \\
        -1 ${reads1} \\
        -2 ${reads2} \\
        -n ${task.cpus} \\
        -o ./${prefix}/${prefix} \\
        \$DEFAULT_ARGS \\
        ${args} || {
        echo "ERROR: RGI bwt failed" >&2
        exit 1
    }

    # Remove intermediate files if they exist
    find "./${prefix}" -type f \\( -name "*.bam" -o -name "*.bam.bai" \\) -delete

    # Verify output files were created
    if [ ! -d "${prefix}" ] || [ -z "\$(ls -A ${prefix} 2>/dev/null)" ]; then
        echo "ERROR: RGI bwt did not produce expected output files" >&2
        exit 1
    fi

    # Clean up symlink
    [ -L "./localDB" ] && rm -f "./localDB"

    # Get tool version with error handling
    RGI_VERSION=\$(rgi main --version 2>/dev/null || echo "unknown")

    # Create versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$RGI_VERSION)
        rgi-database: \$(echo \$DB_VERSION)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.json
    touch ${prefix}/${prefix}.txt
    touch ${prefix}/${prefix}.allele_mapping_data.txt
    touch ${prefix}/${prefix}.gene_mapping_data.txt
    touch ${prefix}/${prefix}.artifacts_mapping_stats.txt
    touch ${prefix}/${prefix}.overall_mapping_stats.txt
    touch ${prefix}/${prefix}.reference_mapping_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: "rgi-stub"
        rgi-database: "4.0.1"
    END_VERSIONS
    """
}