process MERGE_READBASED_POSTPROCESS_INPUTS {
    tag "Aggregate read-based post-processing inputs: ${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.12 conda-forge::pandas=2.2.3"

    input:
    tuple val(meta), path(merged_inputs, stageAs: 'merged_inputs/*', arity: '1..*')

    output:
    tuple val(meta), path('merged_sourmash_yacht.csv'), emit: merged_csv
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail

    declare -a input_files=()
    mapfile -d '' -t input_files < <(
        find -L merged_inputs -mindepth 1 -maxdepth 1 -type f -name '*.csv' -print0 | sort -z
    )
    if [[ "\${#input_files[@]}" -eq 0 ]]; then
        echo 'ERROR: No single-sketch merged CSV files were staged for aggregate post-processing.' >&2
        exit 1
    fi

    python - merged_sourmash_yacht.csv "\${input_files[@]}" <<'PY'
import sys
from pathlib import Path

import pandas as pd

target = Path(sys.argv[1])
input_paths = [Path(path) for path in sys.argv[2:]]
frames = []
for input_path in input_paths:
    try:
        frame = pd.read_csv(input_path, dtype=object)
    except pd.errors.EmptyDataError as exc:
        raise SystemExit(f'ERROR: Aggregate input is empty: {input_path}') from exc
    if frame.columns.empty:
        raise SystemExit(f'ERROR: Aggregate input has no header: {input_path}')
    if 'WGS_ID' not in frame.columns:
        raise SystemExit(f'ERROR: Aggregate input lacks required WGS_ID column: {input_path}')
    if frame['WGS_ID'].isna().any() or frame['WGS_ID'].astype(str).str.strip().eq('').any():
        raise SystemExit(f'ERROR: Aggregate input contains blank WGS_ID values: {input_path}')
    frames.append(frame)

pd.concat(frames, ignore_index=True, sort=False).to_csv(target, index=False)
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch merged_sourmash_yacht.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub_version"
        pandas: "stub_version"
    END_VERSIONS
    """
}
