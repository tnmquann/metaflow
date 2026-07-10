process POSTPROCESS_READBASED {
    tag "${meta.id}:${mode}"
    label 'process_low'

    conda "conda-forge::python=3.12 conda-forge::numpy=2.2.5 conda-forge::pandas=2.2.3"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"

    input:
    tuple val(meta), val(output_prefix), val(mode), env('POSTPROCESS_OPTIONS'), path(merged_sourmash_yacht, stageAs: 'merged_sourmash_yacht.csv'), path(rgi_results, stageAs: 'rgi_inputs/*', arity: '0..*')
    path postprocess_script, stageAs: 'post_processing.py'

    output:
    tuple val(meta), path("${output_prefix}_phyloseq"), emit: phyloseq, optional: true
    tuple val(meta), path("${output_prefix}_taxburst"), emit: taxburst, optional: true
    tuple val(meta), path("${output_prefix}_rgi_bwt"), emit: rgi_bwt, optional: true
    path "versions.yml", emit: versions

    script:
    def run_phyloseq = mode == 'phyloseq' || mode == 'all'
    def run_taxburst = mode == 'taxburst' || mode == 'all'
    def run_rgi_bwt = (mode == 'rgi_bwt' || mode == 'all') && rgi_results
    """
    set -euo pipefail

    mkdir -p postprocess_tmp

    run_postprocess() {
        local output_format="\$1"
        POSTPROCESS_OUTPUT_FORMAT="\${output_format}" python - <<'PY'
    import os
    import shlex
    import subprocess
    import sys

    extra_args = shlex.split(os.environ.get("POSTPROCESS_OPTIONS", ""))
    reserved = ("--input", "--input_format", "--output", "--output_format")
    for argument in extra_args:
        if argument in reserved or any(argument.startswith(f"{flag}=") for flag in reserved):
            raise SystemExit(
                f"postprocess_options cannot override the module-managed argument: {argument}"
            )

    command = [
        sys.executable,
        "post_processing.py",
        *extra_args,
        "--input",
        "merged_sourmash_yacht.csv",
        "--input_format",
        "raw_metaflow",
        "--output",
        "postprocess_tmp",
        "--output_format",
        os.environ["POSTPROCESS_OUTPUT_FORMAT"],
    ]
    subprocess.run(command, check=True)
    PY
    }

    if ${run_phyloseq}; then
        run_postprocess phyloseq
        test -d postprocess_tmp/merged_sourmash_yacht.csv.phyloseq
        mv postprocess_tmp/merged_sourmash_yacht.csv.phyloseq "${output_prefix}_phyloseq"
    fi

    if ${run_taxburst}; then
        run_postprocess taxburst
        test -d postprocess_tmp/merged_sourmash_yacht.csv.taxburst
        mv postprocess_tmp/merged_sourmash_yacht.csv.taxburst "${output_prefix}_taxburst"
    fi

    if ${run_rgi_bwt}; then
        if [[ ! -d rgi_inputs ]] || ! find rgi_inputs -mindepth 1 -print -quit | grep -q .; then
            echo "ERROR: RGI-BWT post-processing was requested, but no staged RGI-BWT results were provided." >&2
            exit 1
        fi
        mkdir -p "${output_prefix}_rgi_bwt"
        cp -R rgi_inputs/. "${output_prefix}_rgi_bwt/"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    def run_phyloseq = mode == 'phyloseq' || mode == 'all'
    def run_taxburst = mode == 'taxburst' || mode == 'all'
    def run_rgi_bwt = (mode == 'rgi_bwt' || mode == 'all') && rgi_results
    """
    if ${run_phyloseq}; then
        mkdir -p "${output_prefix}_phyloseq"
        touch "${output_prefix}_phyloseq/tax_table.csv"
        touch "${output_prefix}_phyloseq/otu_table.csv"
        touch "${output_prefix}_phyloseq/sample_data.csv"
    fi

    if ${run_taxburst}; then
        mkdir -p "${output_prefix}_taxburst"
        touch "${output_prefix}_taxburst/merged_sourmash_yacht.csv.taxburst.csv"
        touch "${output_prefix}_taxburst/sample.taxburst.split.csv"
    fi

    if ${run_rgi_bwt}; then
        mkdir -p "${output_prefix}_rgi_bwt"
        if [[ -d rgi_inputs ]] && find rgi_inputs -mindepth 1 -print -quit | grep -q .; then
            cp -R rgi_inputs/. "${output_prefix}_rgi_bwt/"
        else
            touch "${output_prefix}_rgi_bwt/stub.gene_mapping_data.txt"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub_version"
        numpy: "stub_version"
        pandas: "stub_version"
    END_VERSIONS
    """
}
