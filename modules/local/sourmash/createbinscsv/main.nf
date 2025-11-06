process CREATE_BINS_CSV {
    tag "Creating CSV for bins"
    label 'process_low'

    publishDir "${params.bins_csv_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(bins_dir)

    output:
    tuple val(meta), path("*.csv"), emit: bins_csv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}" ?: "bins"
    """
    #!/bin/bash
    
    # Create CSV header
    echo "name,genome_filename,protein_filename" > ${prefix}_mag.csv
    
    # Process all bin files directly in the working directory
    # Nextflow stages files as symlinks in the current directory or subdirectories
    for bin_file in *.fa *.fasta *.fa.gz *.fasta.gz; do
        # Check if file exists (glob might not match anything)
        if [ -f "\${bin_file}" ]; then
            # Extract just the filename without path and extensions (including .gz)
            bin_name=\$(basename "\${bin_file}" | sed 's/\\.fa\\.gz\$//' | sed 's/\\.fasta\\.gz\$//' | sed 's/\\.fa\$//' | sed 's/\\.fasta\$//')
            # Get full path to the bin file (resolve symlinks)
            full_path=\$(realpath "\${bin_file}")
            # Add entry to CSV (protein_filename column left empty as per requirement)
            echo "\${bin_name},\${full_path}," >> ${prefix}_mag.csv
        fi
    done
    
    # Also check in subdirectories if any
    if [ -d "${bins_dir}" ]; then
        find "${bins_dir}" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fasta.gz" \\) | while read bin_file; do
            # Extract just the filename without path and extensions (including .gz)
            bin_name=\$(basename "\${bin_file}" | sed 's/\\.fa\\.gz\$//' | sed 's/\\.fasta\\.gz\$//' | sed 's/\\.fa\$//' | sed 's/\\.fasta\$//')
            # Get full path to the bin file
            full_path=\$(realpath "\${bin_file}")
            # Add entry to CSV (protein_filename column left empty as per requirement)
            echo "\${bin_name},\${full_path}," >> ${prefix}_mag.csv
        done
    fi

    # Check if any bins were found
    if [ \$(wc -l < ${prefix}_mag.csv) -eq 1 ]; then
        echo "ERROR: No .fa, .fasta, .fa.gz, or .fasta.gz files found"
        echo "Working directory contents:"
        ls -lah
        exit 1
    else
        echo "Found \$((\$(wc -l < ${prefix}_mag.csv) - 1)) bin files"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n 1 | cut -d ' ' -f 4)
        find: \$(find --version 2>&1 | head -n 1 | cut -d ' ' -f 4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.assembler}-${meta.binner ?: 'binrefine'}-${meta.binrefine ?: 'rawbins'}-${meta.id}"
    """
    echo "name,genome_filename,protein_filename" > ${prefix}_mag.csv
    echo "test_bin,/path/to/test.fa," >> ${prefix}_mag.csv
    echo "test_bin_gz,/path/to/test.fa.gz," >> ${prefix}_mag.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub_version"
        find: "stub_version"
    END_VERSIONS
    """
}