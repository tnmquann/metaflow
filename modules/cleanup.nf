process CLEANUP {
    tag "Cleaning temporary files"
    
    input:
    val cleanup_trigger
    
    script:
    """
    #!/bin/bash
    set -e  # Exit on error

    if [ "${params.cleanup}" = "true" ]; then
        echo "Cleaning up temporary files and directories..."
        
        # Create a list of directories to clean
        cleanup_dirs=(
            "${workflow.workDir}"
            "${params.base_dir}/intermediate_files"
            "${params.base_dir}/work"
        )

        # Clean each directory if it exists
        for dir in "\${cleanup_dirs[@]}"; do
            if [ -d "\$dir" ]; then
                echo "Removing directory: \$dir"
                rm -rf "\$dir" || {
                    echo "Warning: Failed to remove \$dir"
                    # Continue even if one directory fails
                    continue
                }
            fi
        done

        # Clean temporary files
        find . -maxdepth 1 -type f -name ".command.*" -delete 2>/dev/null || true
        find . -maxdepth 1 -type f -name ".nextflow.*" -delete 2>/dev/null || true
        
        echo "Directory cleanup completed"
    else
        echo "Cleanup is disabled. Skipping..."
    fi
    """
}
