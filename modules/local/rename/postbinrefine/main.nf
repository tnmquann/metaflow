process RENAME_POSTBINREFINE {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.assembler}-*Refined-*.fa", includeInputs: true), optional: true, emit: refined_bins
    tuple val(meta), path("${meta.assembler}-${meta.binrefine}Unbinned-${meta.id}.fa"), optional: true, emit: refined_unbins
    tuple val(meta), path("${meta.assembler}-${meta.id}-${meta.binrefine}-final_bins_quality_reports.tsv"), optional: true, emit: qc_reports
    tuple val(meta), path("${meta.assembler}-${meta.id}-${meta.binrefine}.log"), optional: true, emit: log_file
    path "versions.yml", emit: versions

    script:
    """
    # Handle quality reports
    if [[ -f final_bins_quality_reports.tsv ]]; then
        mv final_bins_quality_reports.tsv "${meta.assembler}-${meta.id}-${meta.binrefine}-final_bins_quality_reports.tsv"
    fi

    # Handle log file
    if [[ -f "${meta.id}.log" ]]; then
        mv "${meta.id}.log" "${meta.assembler}-${meta.id}-${meta.binrefine}.log"
    fi


    # Handle unbinned contigs
    if [[ -f unbinned.fa ]]; then
        mv unbinned.fa "${meta.assembler}-${meta.binrefine}Unbinned-${meta.id}.fa"
    fi

    if [[ -f bin_*.fa ]]; then
        for binfile in "bin_*.fa"; do
            mv "\$binfile" "${meta.assembler}-${meta.binrefine}-${meta.id}-Refined-\$binfile"
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(echo \$(mv --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.assembler}-*Refined-${meta.id}-bin_1.fa"
    touch "${meta.assembler}-${meta.binrefine}Unbinned-${meta.id}.fa"
    touch "${meta.assembler}-${meta.id}-${meta.binrefine}-final_bins_quality_reports.tsv"
    touch "${meta.assembler}-${meta.id}-${meta.binrefine}.log"
    touch "versions.yml"
    """
}