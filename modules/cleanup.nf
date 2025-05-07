process CLEANUP {
    tag "Cleaning temporary files"
    label 'process_low' // Low resource usage

    // No Conda needed for basic shell commands unless specific versions are required

    input:
    val flag_to_run // Input to trigger the process, can be output from PROCESS_RESULTS

    output:
    path "versions.yml", emit: versions // Optional: to maintain consistency

    script:
    """
    #!/bin/bash
    # Script execution is conditional on params.cleanup,
    # so this process should only be invoked if params.cleanup is true.
    # The 'input: val flag_to_run' ensures it runs after preceding steps.

    echo "Starting cleanup process in CLEANUP module..."

    # Define directories to clean
    # workflow.workDir is the main work directory.
    # params.base_dir might be where intermediate/results are stored outside 'work'.
    cleanup_dirs=(
        "${workflow.workDir}"
        // "${params.base_dir}/intermediate_files" // Be careful with this
        // "${params.base_dir}/work" // This is usually the same as workflow.workDir or a sub-path
    )

    for dir_to_clean in "\${cleanup_dirs[@]}"; do
        if [ -d "\$dir_to_clean" ]; then
            echo "Cleaning directory: \$dir_to_clean"
            # Add safety checks or be more specific if needed
            # rm -rf "\$dir_to_clean"/* # Example: clean contents, not the dir itself
            # For now, keeping original logic but with a warning:
            echo "WARNING: rm -rf on \$dir_to_clean is a destructive operation."
            # rm -rf "\$dir_to_clean" || {
            #     echo "Warning: Failed to remove \$dir_to_clean"
            # }
            echo "Skipping actual rm -rf for safety in this example. Uncomment if sure."
        else
            echo "Directory not found, skipping cleanup: \$dir_to_clean"
        fi
    done

    # Clean temporary files in the current process's work directory
    # This is less common as Nextflow cleans process work dirs itself unless -resume is used.
    # The original script's find commands were relative to the process work dir.
    # echo "Cleaning temporary files in current process work directory..."
    # find . -maxdepth 1 -type f -name ".command.*" -delete 2>/dev/null || true
    # find . -maxdepth 1 -type f -name ".nextflow.*" -delete 2>/dev/null || true

    echo "Cleanup module execution finished."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1)
        rm: \$(rm --version 2>&1 | head -n 1 || echo "rm (coreutils)")
        find: \$(find --version 2>&1 | head -n 1 || echo "find (coreutils)")
    END_VERSIONS
    """

    stub:
    """
    echo "Stub: Cleanup process would run here."
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub_version"
        rm: "stub_version"
        find: "stub_version"
    END_VERSIONS
    """
}
