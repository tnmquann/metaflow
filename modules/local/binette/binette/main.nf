process BINETTE_BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/binette:1.2.0--pyhdfd78af_0':
        'quay.io/biocontainers/binette:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(bins), path(proteins)
    path(checkm2_db)

    output:
    tuple val(meta), path("${prefix}_final_bins_quality_reports.tsv")              , emit: quality_reports
    tuple val(meta), path("${prefix}_final_bins/")                 , optional: true , emit: bins
    tuple val(meta), path("${prefix}_final_contig_to_bin.tsv")                     , emit: contig2bin
    tuple val(meta), path("${prefix}_input_bins_quality_reports/")                 , emit: input_reports
    tuple val(meta), path("${prefix}.log")                                         , emit: log
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_arg = proteins ? "--proteins $proteins" : ""
    def checkm2_db_arg = checkm2_db ? "--checkm2_db $checkm2_db" : ""
    
    // Handle two different input modes like DAS_Tool
    def bins_arg = ""
    if (meta.input_mode == "contig2bin_tables") {
        // bins is a list of TSV files for contig2bin tables
        def bin_list = bins instanceof List ? bins.join(" ") : "$bins"
        bins_arg = "--contig2bin_tables $bin_list"
    } else if (meta.input_mode == "bin_dirs") {
        // bins is a list of directories containing bin FASTA files  
        def bin_list = bins instanceof List ? bins.join(" ") : "$bins"
        bins_arg = "--bin_dirs $bin_list"
    } else {
        error "meta.input_mode must be either 'contig2bin_tables' or 'bin_dirs'"
    }

    """
    binette \\
        $bins_arg \\
        --contigs $contigs \\
        $proteins_arg \\
        $checkm2_db_arg \\
        --threads $task.cpus \\
        --outdir ${prefix}_binette_output \\
        --prefix $prefix \\
        $args \\
        2>&1 | tee ${prefix}.log

    # Move outputs to expected locations
    mv ${prefix}_binette_output/final_bins_quality_reports.tsv ${prefix}_final_bins_quality_reports.tsv
    mv ${prefix}_binette_output/final_contig_to_bin.tsv ${prefix}_final_contig_to_bin.tsv
    mv ${prefix}_binette_output/input_bins_quality_reports ${prefix}_input_bins_quality_reports
    
    # Only move final_bins directory if it exists (depends on --write-fasta-bins setting)
    if [ -d "${prefix}_binette_output/final_bins" ]; then
        mv ${prefix}_binette_output/final_bins ${prefix}_final_bins
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(binette --version 2>&1 | sed 's/Binette //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_final_bins_quality_reports.tsv
    mkdir ${prefix}_final_bins
    touch ${prefix}_final_bins/bin1.fa
    touch ${prefix}_final_contig_to_bin.tsv
    mkdir ${prefix}_input_bins_quality_reports
    touch ${prefix}_input_bins_quality_reports/report.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(binette --version 2>&1 | sed 's/Binette //' || echo "1.2.0")
    END_VERSIONS
    """
}