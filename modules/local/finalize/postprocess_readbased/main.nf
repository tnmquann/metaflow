process POSTPROCESS_READBASED {
    tag "${meta.id}:${mode}"
    label 'process_low'

    conda "conda-forge::python=3.12 conda-forge::numpy=2.2.5 conda-forge::pandas=2.2.3"
    container "community.wave.seqera.io/library/sourmash_pandas_sourmash_plugin_branchwater:477d25f3865da957"

    input:
    tuple val(meta), val(output_prefix), val(mode), val(single_sample_layout), val(rgi_gene_mapping_only), env('POSTPROCESS_OPTIONS'), path(merged_sourmash_yacht, stageAs: 'merged_sourmash_yacht.csv'), path(rgi_results, stageAs: 'rgi_inputs/*', arity: '0..*')
    path postprocess_script, stageAs: 'post_processing.py'

    output:
    tuple val(meta), path("${output_prefix}_default"), emit: default_output, optional: true
    tuple val(meta), path("${output_prefix}_phyloseq"), emit: phyloseq, optional: true
    tuple val(meta), path("${output_prefix}_taxburst"), emit: taxburst, optional: true
    tuple val(meta), path("${output_prefix}_rgi_bwt"), emit: rgi_bwt, optional: true
    path "versions.yml", emit: versions

    script:
    def run_default = mode == 'all'
    def run_phyloseq = mode == 'phyloseq' || mode == 'all'
    def run_taxburst = mode == 'taxburst' || mode == 'all'
    def run_rgi_bwt = (mode == 'rgi_bwt' || mode == 'all') && rgi_results
    """
    set -euo pipefail

    run_postprocess_command() {
        local output_format="\$1"
        local mode_tmp="\$2"

        POSTPROCESS_OUTPUT_FORMAT="\${output_format}" POSTPROCESS_TMP_DIR="\${mode_tmp}" python - <<'PY'
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
        os.environ["POSTPROCESS_TMP_DIR"],
        "--output_format",
        os.environ["POSTPROCESS_OUTPUT_FORMAT"],
    ]
    subprocess.run(command, check=True)
    PY
    }

    run_postprocess() {
        local output_format="\$1"
        local final_output="\$2"
        local mode_tmp="postprocess_tmp/\${output_format}"
        local -a generated_dirs=()

        rm -rf "\${mode_tmp}"
        mkdir -p "\${mode_tmp}"

        run_postprocess_command "\${output_format}" "\${mode_tmp}"

        mapfile -d '' -t generated_dirs < <(
            find "\${mode_tmp}" -mindepth 1 -maxdepth 1 -type d -name "*.\${output_format}" -print0
        )
        if [[ "\${#generated_dirs[@]}" -ne 1 ]]; then
            echo "ERROR: Expected exactly one \${output_format} output directory in \${mode_tmp}; found \${#generated_dirs[@]}." >&2
            find "\${mode_tmp}" -maxdepth 2 -print >&2
            exit 1
        fi

        mv "\${generated_dirs[0]}" "\${final_output}"
    }

    run_default_postprocess() {
        local final_output="\$1"
        local mode_tmp='postprocess_tmp/default'
        local -a generated_files=()
        local -a remaining_outputs=()
        local default_csv_name="${output_prefix}.processed_metaflow.csv"

        rm -rf "\${mode_tmp}"
        mkdir -p "\${mode_tmp}"

        run_postprocess_command default "\${mode_tmp}"

        mapfile -d '' -t generated_files < <(
            find "\${mode_tmp}" -mindepth 1 -maxdepth 1 -type f -name '*.processed_metaflow.csv' -print0
        )
        if [[ "\${#generated_files[@]}" -ne 1 ]]; then
            echo "ERROR: Expected exactly one default processed CSV in \${mode_tmp}; found \${#generated_files[@]}." >&2
            find "\${mode_tmp}" -maxdepth 2 -print >&2
            exit 1
        fi

        mkdir -p "\${final_output}"
        if [[ -e "\${final_output}/\${default_csv_name}" ]]; then
            echo "ERROR: Default processed CSV destination already exists: \${final_output}/\${default_csv_name}." >&2
            exit 1
        fi
        mv "\${generated_files[0]}" "\${final_output}/\${default_csv_name}"
        shopt -s dotglob nullglob
        remaining_outputs=("\${mode_tmp}"/*)
        if [[ "\${#remaining_outputs[@]}" -gt 0 ]]; then
            mv "\${remaining_outputs[@]}" "\${final_output}/"
        fi
    }

    if ${run_default}; then
        run_default_postprocess "${output_prefix}_default"
    fi

    if ${run_phyloseq}; then
        run_postprocess phyloseq "${output_prefix}_phyloseq"
    fi

    if ${run_taxburst}; then
        run_postprocess taxburst "${output_prefix}_taxburst"
    fi

    if ${run_rgi_bwt}; then
        if [[ ! -d rgi_inputs ]] || ! find -L rgi_inputs -mindepth 1 -print -quit | grep -q .; then
            echo "ERROR: RGI-BWT post-processing was requested, but no staged RGI-BWT results were provided." >&2
            exit 1
        fi
        mkdir -p "${output_prefix}_rgi_bwt"
        if ${rgi_gene_mapping_only}; then
            declare -a gene_mapping_files=()
            mapfile -d '' -t gene_mapping_files < <(
                find -L rgi_inputs -type f -name '*.gene_mapping_data.txt' -print0 | sort -z
            )
            if [[ "\${#gene_mapping_files[@]}" -eq 0 ]]; then
                echo 'ERROR: RGI-BWT post-processing expected at least one *.gene_mapping_data.txt file.' >&2
                exit 1
            fi
            for gene_mapping_file in "\${gene_mapping_files[@]}"; do
                target_file="${output_prefix}_rgi_bwt/\$(basename "\${gene_mapping_file}")"
                if [[ -e "\${target_file}" ]]; then
                    echo "ERROR: Multiple RGI-BWT gene mapping files have the same filename: \$(basename "\${gene_mapping_file}")." >&2
                    exit 1
                fi
                cp "\${gene_mapping_file}" "\${target_file}"
            done
        else
            cp -R rgi_inputs/. "${output_prefix}_rgi_bwt/"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    def run_default = mode == 'all'
    def run_phyloseq = mode == 'phyloseq' || mode == 'all'
    def run_taxburst = mode == 'taxburst' || mode == 'all'
    def run_rgi_bwt = (mode == 'rgi_bwt' || mode == 'all') && rgi_results
    """
    if ${run_default}; then
        mkdir -p "${output_prefix}_default"
        touch "${output_prefix}_default/${output_prefix}.processed_metaflow.csv"
    fi

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
        if [[ -d rgi_inputs ]] && find -L rgi_inputs -mindepth 1 -print -quit | grep -q .; then
            if ${rgi_gene_mapping_only}; then
                declare -a gene_mapping_files=()
                mapfile -d '' -t gene_mapping_files < <(
                    find -L rgi_inputs -type f -name '*.gene_mapping_data.txt' -print0 | sort -z
                )
                if [[ "\${#gene_mapping_files[@]}" -eq 0 ]]; then
                    echo 'ERROR: RGI-BWT post-processing expected at least one *.gene_mapping_data.txt file.' >&2
                    exit 1
                fi
                for gene_mapping_file in "\${gene_mapping_files[@]}"; do
                    target_file="${output_prefix}_rgi_bwt/\$(basename "\${gene_mapping_file}")"
                    if [[ -e "\${target_file}" ]]; then
                        echo "ERROR: Multiple RGI-BWT gene mapping files have the same filename: \$(basename "\${gene_mapping_file}")." >&2
                        exit 1
                    fi
                    cp "\${gene_mapping_file}" "\${target_file}"
                done
            else
                cp -R rgi_inputs/. "${output_prefix}_rgi_bwt/"
            fi
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
