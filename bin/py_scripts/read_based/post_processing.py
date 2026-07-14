#!/usr/bin/env python3
"""Process plain-text metagenomics tabular files."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import logging
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


LOGGER = logging.getLogger("metagenomic_processor")


GATHER_REQUIRED_COLUMNS = [
    "query_name",
    "scaled",
    "ksize",
]

RAW_METAFLOW_STANDARD_COLUMNS = [
    "sample_name",
    "lineage",
    "query_md5",
    "query_filename",
    "ksize",
    "scaled",
    "match_name",
    "unweighted_fraction",
    "f_weighted_at_rank",
    "bp_match_at_rank",
    "total_weighted_hashes",
    "n_unique_kmers",
    "p_vals",
    "actual_confidence_with_coverage",
    "alt_confidence_mut_rate_with_coverage",
    "min_coverage",
    "source_file",
]

RAW_METAFLOW_TAXBURST_EXTRA_COLUMNS = [
    "n_unique_weighted_found",
    "sum_weighted_found",
    "median_abund",
]

RAW_METAFLOW_NUMERIC_COLUMNS = [
    "f_unique_to_query",
    "average_abund",
    "median_abund",
    "unique_intersect_bp",
    "query_bp",
    "scaled",
    "relative_abundance",
    "total_weighted_hashes",
    "n_unique_weighted_found",
    "sum_weighted_found",
]

RAW_METAFLOW_STRICT_GROUP_METADATA_COLUMNS = [
    "query_bp",
    "total_weighted_hashes",
]

RAW_METAFLOW_SOFT_GROUP_METADATA_COLUMNS = [
    "p_vals",
    "actual_confidence_with_coverage",
    "alt_confidence_mut_rate_with_coverage",
]

RAW_METAFLOW_REFERENCED_COLUMNS = [
    "WGS_ID",
    "lineage",
    "query_md5",
    "query_filename",
    "ksize",
    "scaled",
    "organism_name",
    "f_unique_to_query",
    "unique_intersect_bp",
    "query_bp",
    "total_weighted_hashes",
    "average_abund",
    "median_abund",
    "relative_abundance",
    "p_vals",
    "actual_confidence_with_coverage",
    "alt_confidence_mut_rate_with_coverage",
    "file_name",
    "n_unique_weighted_found",
    "sum_weighted_found",
]

RAW_METAFLOW_PREFERRED_GROUP_COLUMNS = [
    "__sample_name",
    "file_name",
    "query_md5",
    "query_filename",
    "ksize",
    "scaled",
]

RAW_METAFLOW_MINIMUM_INFORMATIVE_COLUMNS = [
    "organism_name",
    "WGS_ID",
    "query_md5",
    "query_filename",
    "relative_abundance",
    "file_name",
]

SUPPORTED_EXTENSIONS = {".csv", ".tsv", ".txt"}
INVALID_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
MIN_COVERAGE_PATTERN = re.compile(r"min_coverage([0-9eE.+-]+)")


class ProcessingError(RuntimeError):
    """Raised when input validation or processing fails."""


@dataclass
class OutputTarget:
    mode: str
    path: Path | None


@dataclass
class ProcessedResult:
    dataframe: pd.DataFrame
    source_path: Path
    wgs_ids_complete: bool


@dataclass(frozen=True)
class SplitTarget:
    directory: Path
    filename_prefix: str = ""


@dataclass
class RawMetaflowWarningCollector:
    input_path: Path
    messages: list[str] = field(default_factory=list)
    seen: set[str] = field(default_factory=set)

    def add(self, message: str) -> None:
        if message in self.seen:
            return
        self.seen.add(message)
        self.messages.append(message)

    def emit(self) -> None:
        if not self.messages:
            return
        details = "\n".join(f"  - {message}" for message in self.messages)
        LOGGER.warning("raw_metaflow warnings for %s:\n%s", self.input_path, details)


def emit_warning(
    message: str,
    *,
    collector: RawMetaflowWarningCollector | None = None,
) -> None:
    if collector is None:
        LOGGER.warning("%s", message)
        return
    collector.add(message)


def build_parser() -> argparse.ArgumentParser:
    examples = """Examples:
  python process_metagenomic_microgen2026.py --input "data/*.csv" --input_format gather
  python process_metagenomic_microgen2026.py --input "data/{a,b}.tsv" --input_format raw_metaflow --output results.csv
  python process_metagenomic_microgen2026.py --input merged.csv --input_format raw_metaflow --remove_unclassified --split
  python process_metagenomic_microgen2026.py --input merged.csv --input_format raw_metaflow --min_coverage 0.01 --output_format phyloseq
  python process_metagenomic_microgen2026.py --input merged.csv --input_format raw_metaflow --output_format taxburst
"""
    parser = argparse.ArgumentParser(
        description=(
            "Process plain-text CSV/TSV/TXT metagenomics tabular files for gather "
            "splitting or raw_metaflow standardization."
        ),
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        required=True,
        nargs="+",
        help="One or more input path expressions. Supports shell-style globbing and brace expansion.",
    )
    parser.add_argument(
        "--input_format",
        required=True,
        choices=["gather", "raw_metaflow"],
        help="Declared schema for all input files.",
    )
    parser.add_argument(
        "--output",
        help=(
            "Optional output filename or existing directory. For gather splitting and "
            "phyloseq/taxburst exports, use an existing directory or omit this option."
        ),
    )
    parser.add_argument(
        "--split",
        action="store_true",
        help="For raw_metaflow only: also create one processed CSV per sample_name. taxburst always emits split files.",
    )
    parser.add_argument(
        "--remove_unclassified",
        action="store_true",
        help="For raw_metaflow only: omit synthetic unclassified rows and renormalize f_weighted_at_rank.",
    )
    parser.add_argument(
        "--min_coverage",
        type=float,
        help="For raw_metaflow only: keep rows whose file_name contains min_coverage{value}.",
    )
    parser.add_argument(
        "--use_average",
        action="store_true",
        help="Use average_abund instead of median_abund in abundance formulas.",
    )
    parser.add_argument(
        "--output_format",
        choices=["default", "phyloseq", "taxburst"],
        default="default",
        help="Output layout for raw_metaflow processing. Supported: default, phyloseq, taxburst. Default: default.",
    )
    return parser


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        validate_cli_arguments(args, parser)
        input_paths = expand_input_paths(args.input)
        output_target = determine_output_target(args.output)

        if args.input_format == "gather":
            process_gather_inputs(input_paths, output_target)
        else:
            process_raw_metaflow_inputs(input_paths, output_target, args)
    except ProcessingError as exc:
        LOGGER.error(str(exc))
        return 1
    except KeyboardInterrupt:
        LOGGER.error("Interrupted by user.")
        return 130

    return 0


def validate_cli_arguments(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if args.split and args.input_format != "raw_metaflow":
        parser.error("--split is only valid with --input_format raw_metaflow.")
    if args.remove_unclassified and args.input_format != "raw_metaflow":
        parser.error("--remove_unclassified is only valid with --input_format raw_metaflow.")
    if args.min_coverage is not None and args.input_format != "raw_metaflow":
        parser.error("--min_coverage is only valid with --input_format raw_metaflow.")
    if args.output_format != "default" and args.input_format != "raw_metaflow":
        parser.error(f"--output_format {args.output_format} is only valid with --input_format raw_metaflow.")
    if args.output_format == "phyloseq" and args.split:
        parser.error("--split is not supported with --output_format phyloseq.")
    if args.output_format == "taxburst" and args.remove_unclassified:
        parser.error(
            "--remove_unclassified cannot be combined with --output_format taxburst because taxburst already omits synthetic unclassified rows without renormalizing."
        )
    if args.min_coverage is not None and args.min_coverage < 0:
        parser.error("--min_coverage must be non-negative.")


def expand_input_paths(expressions: Sequence[str]) -> list[Path]:
    expanded: list[Path] = []
    seen: set[Path] = set()

    for expression in expressions:
        matched_for_expression = False
        for brace_expanded in expand_braces(expression):
            literal = os_safe_expanduser(brace_expanded)
            matches = glob.glob(literal)
            if not matches and Path(literal).exists():
                matches = [literal]
            for match in matches:
                path = Path(match).resolve()
                matched_for_expression = True
                if path in seen:
                    continue
                if not path.is_file():
                    raise ProcessingError(f"Matched path is not a file: {path}")
                if path.suffix.lower() not in SUPPORTED_EXTENSIONS:
                    raise ProcessingError(
                        f"Unsupported input extension for {path}. Supported extensions: "
                        f"{', '.join(sorted(SUPPORTED_EXTENSIONS))}."
                    )
                seen.add(path)
                expanded.append(path)
        if not matched_for_expression:
            raise ProcessingError(f"No files matched input expression: {expression}")

    if not expanded:
        raise ProcessingError("No input files matched the provided expressions.")
    return expanded


def os_safe_expanduser(value: str) -> str:
    return str(Path(value).expanduser())


def expand_braces(pattern: str) -> list[str]:
    start = pattern.find("{")
    if start == -1:
        return [pattern]

    depth = 0
    end = -1
    for index in range(start, len(pattern)):
        char = pattern[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                end = index
                break

    if end == -1:
        return [pattern]

    prefix = pattern[:start]
    suffix = pattern[end + 1 :]
    body = pattern[start + 1 : end]
    options = split_brace_options(body)

    expanded: list[str] = []
    for option in options:
        expanded.extend(expand_braces(prefix + option + suffix))
    return expanded


def split_brace_options(body: str) -> list[str]:
    options: list[str] = []
    depth = 0
    current: list[str] = []

    for char in body:
        if char == "," and depth == 0:
            options.append("".join(current))
            current = []
            continue
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
        current.append(char)

    options.append("".join(current))
    return options


def determine_output_target(output_arg: str | None) -> OutputTarget:
    if output_arg is None:
        return OutputTarget(mode="auto", path=None)

    output_path = Path(output_arg).expanduser()
    if output_path.exists() and output_path.is_dir():
        return OutputTarget(mode="directory", path=output_path.resolve())

    parent = output_path.parent if output_path.parent != Path("") else Path(".")
    if not parent.exists():
        raise ProcessingError(f"Output parent directory does not exist: {parent}")

    return OutputTarget(mode="file", path=output_path.resolve())


def detect_delimiter(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ","
    if suffix == ".tsv":
        return "\t"

    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(8192)
    except UnicodeDecodeError as exc:
        raise ProcessingError(f"Failed to read {path} as UTF-8 text: {exc}") from exc

    if not sample.strip():
        raise ProcessingError(f"Input file is empty: {path}")

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except csv.Error:
        first_line = sample.splitlines()[0]
        if "\t" in first_line and first_line.count("\t") >= first_line.count(","):
            return "\t"
        return ","


def read_table(path: Path) -> pd.DataFrame:
    delimiter = detect_delimiter(path)
    try:
        dataframe = pd.read_csv(
            path,
            sep=delimiter,
            encoding="utf-8-sig",
            dtype=object,
            keep_default_na=True,
            na_values=["", "NA", "NaN", "null", "NULL"],
        )
    except pd.errors.EmptyDataError as exc:
        raise ProcessingError(f"Input file is empty: {path}") from exc
    except UnicodeDecodeError as exc:
        raise ProcessingError(f"Failed to decode {path} as UTF-8 text: {exc}") from exc
    except Exception as exc:  # pragma: no cover - defensive parsing guard
        raise ProcessingError(f"Failed to read {path}: {exc}") from exc

    if dataframe.columns.empty:
        raise ProcessingError(f"Input file has no header row: {path}")
    if dataframe.empty:
        raise ProcessingError(f"Input file contains a header but no records: {path}")

    dataframe.columns = [str(column).strip() for column in dataframe.columns]
    return dataframe


def validate_required_columns(dataframe: pd.DataFrame, required_columns: Sequence[str], path: Path) -> None:
    missing = [column for column in required_columns if column not in dataframe.columns]
    if missing:
        raise ProcessingError(
            f"{path} is missing required columns: {', '.join(missing)}"
        )


def add_missing_columns(dataframe: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    for column in columns:
        if column not in dataframe.columns:
            dataframe[column] = pd.NA
    return dataframe


def warn_on_missing_optional_columns(
    dataframe: pd.DataFrame,
    optional_columns: Sequence[str],
    path: Path,
    *,
    context: str,
    collector: RawMetaflowWarningCollector | None = None,
) -> None:
    missing = [column for column in optional_columns if column not in dataframe.columns]
    if missing:
        emit_warning(
            f"Missing non-critical {context} columns: {', '.join(missing)}. "
            "Related output fields will be left blank where possible.",
            collector=collector,
        )


def coerce_numeric_columns(
    dataframe: pd.DataFrame,
    columns: Sequence[str],
    path: Path,
    *,
    allow_missing: bool = True,
) -> pd.DataFrame:
    for column in columns:
        if column not in dataframe.columns:
            continue
        original = dataframe[column]
        numeric = pd.to_numeric(original, errors="coerce")

        nonempty = original.notna() & original.astype(str).str.strip().ne("")
        invalid = nonempty & numeric.isna()
        if invalid.any():
            examples = ", ".join(sorted({str(value) for value in original[invalid].head(5)}))
            raise ProcessingError(
                f"{path} contains non-numeric values in column '{column}': {examples}"
            )
        if not allow_missing and numeric.isna().any():
            raise ProcessingError(
                f"{path} contains missing numeric values in required column '{column}'."
            )
        dataframe[column] = numeric
    return dataframe


def require_non_null_values(dataframe: pd.DataFrame, columns: Sequence[str], path: Path) -> None:
    for column in columns:
        missing = dataframe[column].isna()
        if missing.any():
            raise ProcessingError(
                f"{path} contains missing values in required column '{column}' for {missing.sum()} rows."
            )


def validate_group_key_values(dataframe: pd.DataFrame, columns: Sequence[str], path: Path) -> None:
    for column in columns:
        missing = dataframe[column].isna() | dataframe[column].astype(str).str.strip().eq("")
        if missing.any():
            raise ProcessingError(
                f"{path} contains missing values in required grouping column '{column}' for {missing.sum()} rows."
            )


def resolve_raw_metaflow_required_columns(
    *,
    use_average: bool,
    include_synthetic_unclassified: bool,
    min_coverage: float | None,
) -> list[str]:
    formula_column = "average_abund" if use_average else "median_abund"
    required = [
        formula_column,
        "unique_intersect_bp",
        "scaled",
        "total_weighted_hashes",
    ]
    if include_synthetic_unclassified:
        required.extend(["query_bp", "f_unique_to_query"])
    if min_coverage is not None:
        required.extend(["relative_abundance", "file_name"])
    return list(dict.fromkeys(required))


def resolve_raw_metaflow_non_null_columns(
    *,
    use_average: bool,
    include_synthetic_unclassified: bool,
    min_coverage: float | None,
) -> list[str]:
    required = resolve_raw_metaflow_required_columns(
        use_average=use_average,
        include_synthetic_unclassified=include_synthetic_unclassified,
        min_coverage=min_coverage,
    )
    return required


def validate_raw_metaflow_minimum_columns(dataframe: pd.DataFrame, path: Path) -> None:
    present = [column for column in RAW_METAFLOW_MINIMUM_INFORMATIVE_COLUMNS if column in dataframe.columns]
    if not present:
        raise ProcessingError(
            f"{path} does not contain any recognizable raw_metaflow identification/result columns. "
            f"Expected at least one of: {', '.join(RAW_METAFLOW_MINIMUM_INFORMATIVE_COLUMNS)}"
        )


def split_formula_ready_rows(
    dataframe: pd.DataFrame,
    required_columns: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if dataframe.empty:
        return dataframe.copy(), dataframe.copy()

    readiness_mask = pd.Series(True, index=dataframe.index)
    for column in required_columns:
        readiness_mask = readiness_mask & dataframe[column].notna()

    return dataframe[readiness_mask].copy(), dataframe[~readiness_mask].copy()


def resolve_raw_output_columns(output_format: str) -> list[str]:
    columns = list(RAW_METAFLOW_STANDARD_COLUMNS)
    if output_format == "taxburst":
        columns.extend(RAW_METAFLOW_TAXBURST_EXTRA_COLUMNS)
    return columns


def raw_output_includes_synthetic_unclassified(output_format: str, remove_unclassified: bool) -> bool:
    return output_format != "taxburst" and not remove_unclassified


def raw_output_requires_split_files(output_format: str, split_requested: bool) -> bool:
    return split_requested or output_format == "taxburst"


def frames_have_matching_columns(frames: Sequence[pd.DataFrame]) -> bool:
    if not frames:
        return True
    expected = list(frames[0].columns)
    return all(list(frame.columns) == expected for frame in frames[1:])


def build_sample_name_series(dataframe: pd.DataFrame, input_path: Path) -> pd.Series:
    fallback = pd.Series(input_path.stem, index=dataframe.index, dtype="object")
    if "WGS_ID" not in dataframe.columns:
        return fallback

    values = dataframe["WGS_ID"].astype("string")
    nonempty = values.notna() & values.str.strip().ne("")
    return values.where(nonempty, other=input_path.stem).astype("object")


def column_has_meaningful_values(series: pd.Series) -> bool:
    if pd.api.types.is_numeric_dtype(series):
        return series.notna().any()

    text = series.astype("string")
    return bool((text.notna() & text.str.strip().ne("")).any())


def resolve_raw_metaflow_group_columns(
    dataframe: pd.DataFrame,
    input_path: Path,
    *,
    include_file_name: bool = True,
    warn_on_reduced: bool = True,
    collector: RawMetaflowWarningCollector | None = None,
) -> list[str]:
    available: list[str] = []
    reduced: list[str] = []

    for column in RAW_METAFLOW_PREFERRED_GROUP_COLUMNS:
        if column == "__sample_name__":
            available.append(column)
            continue
        if column == "file_name" and not include_file_name:
            continue
        if column in dataframe.columns and column_has_meaningful_values(dataframe[column]):
            available.append(column)
        else:
            reduced.append(column)

    if reduced and warn_on_reduced:
        emit_warning(
            f"raw_metaflow grouping is using a reduced key because these columns are missing or empty: {', '.join(reduced)}",
            collector=collector,
        )

    return available


def build_passthrough_precomputed_rows(dataframe: pd.DataFrame, *, output_format: str) -> pd.DataFrame:
    standardized = pd.DataFrame(
        {
            "sample_name": dataframe["__sample_name"],
            "lineage": dataframe["lineage"],
            "query_md5": dataframe["query_md5"],
            "query_filename": dataframe["query_filename"],
            "ksize": dataframe["ksize"],
            "scaled": dataframe["scaled"],
            "match_name": dataframe["organism_name"],
            "unweighted_fraction": dataframe["f_unique_to_query"],
            "f_weighted_at_rank": dataframe["relative_abundance"],
            "bp_match_at_rank": dataframe["unique_intersect_bp"],
            "total_weighted_hashes": dataframe["total_weighted_hashes"],
            "n_unique_kmers": pd.Series(pd.NA, index=dataframe.index, dtype="object"),
            "p_vals": dataframe["p_vals"],
            "actual_confidence_with_coverage": dataframe["actual_confidence_with_coverage"],
            "alt_confidence_mut_rate_with_coverage": dataframe["alt_confidence_mut_rate_with_coverage"],
            "min_coverage": dataframe["file_name"],
            "source_file": dataframe["__provenance"],
        }
    )
    if output_format == "taxburst":
        standardized["n_unique_weighted_found"] = dataframe["n_unique_weighted_found"]
        standardized["sum_weighted_found"] = dataframe["sum_weighted_found"]
        standardized["median_abund"] = dataframe["median_abund"]
    return standardized.reindex(columns=resolve_raw_output_columns(output_format))


def drop_blank_raw_metaflow_records(
    dataframe: pd.DataFrame,
    input_path: Path,
    *,
    collector: RawMetaflowWarningCollector | None = None,
) -> pd.DataFrame:
    file_name_blank = dataframe["file_name"].isna()
    if "file_name" in dataframe.columns:
        file_name_blank = file_name_blank | dataframe["file_name"].astype("string").str.strip().eq("")
    if file_name_blank.any():
        emit_warning(
            f"{int(file_name_blank.sum())} raw_metaflow rows have blank file_name and will be dropped before downstream processing.",
            collector=collector,
        )
        dataframe = dataframe.loc[~file_name_blank].copy()

    relative_abundance_blank = dataframe["relative_abundance"].isna()
    if relative_abundance_blank.any():
        emit_warning(
            f"{int(relative_abundance_blank.sum())} raw_metaflow rows have blank relative_abundance after blank file_name rows were removed and will be dropped before downstream processing.",
            collector=collector,
        )
        dataframe = dataframe.loc[~relative_abundance_blank].copy()

    return dataframe


def process_gather_inputs(input_paths: Sequence[Path], output_target: OutputTarget) -> None:
    merge_outputs = output_target.mode == "file"
    merge_frames: list[pd.DataFrame] = []
    grouped_outputs: list[tuple[Path, pd.DataFrame]] = []

    for input_path in input_paths:
        dataframe = read_table(input_path)
        validate_required_columns(dataframe, GATHER_REQUIRED_COLUMNS, input_path)
        validate_group_key_values(dataframe, ["query_name", "scaled", "ksize"], input_path)

        groups = list(dataframe.groupby(["query_name", "scaled", "ksize"], dropna=False, sort=False))
        if len(groups) == 1:
            LOGGER.warning(
                "Only one gather group was found in %s; still writing the standard grouped output.",
                input_path,
            )

        if merge_outputs:
            merge_frame = dataframe.copy()
            merge_frame["source_file"] = str(input_path)
            merge_frames.append(merge_frame)
            continue

        for group_key, group_frame in groups:
            query_name, scaled, ksize = group_key
            filename = (
                f"{sanitize_filename(format_value(query_name))}.s"
                f"{sanitize_filename(format_value(scaled))}k"
                f"{sanitize_filename(format_value(ksize))}.gather.csv"
            )
            output_path = resolve_gather_output_path(input_path, output_target, filename)
            grouped_outputs.append((output_path, group_frame.copy()))

    if merge_outputs:
        if not frames_have_matching_columns(merge_frames):
            LOGGER.warning(
                "Gather inputs do not share identical optional columns; merged output will use the union of columns and leave missing values blank."
            )
        merged = pd.concat(merge_frames, ignore_index=True, sort=False)
        if len(input_paths) > 1:
            LOGGER.warning("Multiple input files are being merged into one output CSV.")
        assert output_target.path is not None
        write_csv(merged, output_target.path)
        LOGGER.info("Wrote %s", output_target.path)
        return

    ensure_unique_output_paths([path for path, _ in grouped_outputs])

    for output_path, group_frame in grouped_outputs:
        write_csv(group_frame, output_path)
        LOGGER.info("Wrote %s", output_path)


def resolve_gather_output_path(input_path: Path, output_target: OutputTarget, filename: str) -> Path:
    if output_target.mode == "auto":
        return input_path.parent / filename
    if output_target.mode == "directory":
        assert output_target.path is not None
        return output_target.path / filename
    assert output_target.path is not None
    return output_target.path


def process_raw_metaflow_inputs(
    input_paths: Sequence[Path],
    output_target: OutputTarget,
    args: argparse.Namespace,
) -> None:
    merge_outputs = len(input_paths) > 1 and output_target.mode == "file" and args.output_format == "default"
    split_enabled = raw_output_requires_split_files(args.output_format, args.split)

    if args.output_format in {"phyloseq", "taxburst"} and output_target.mode == "file":
        raise ProcessingError(
            f"--output_format {args.output_format} produces multiple files. Use an existing directory or omit --output."
        )

    processed_results = [
        process_raw_metaflow_file(
            input_path=path,
            use_average=args.use_average,
            remove_unclassified=args.remove_unclassified,
            min_coverage=args.min_coverage,
            output_format=args.output_format,
            suppress_remove_unclassified_warning=args.output_format == "phyloseq",
        )
        for path in input_paths
    ]
    output_labels = build_raw_output_labels([result.source_path for result in processed_results], output_target)
    split_targets = build_split_targets(processed_results, output_target, args.output_format, output_labels)

    if merge_outputs:
        LOGGER.warning("Multiple input files are being merged into one output CSV.")
        merged = merge_processed_raw_results(processed_results)
        assert output_target.path is not None
        write_csv(merged, output_target.path)
        LOGGER.info("Wrote %s", output_target.path)
        if split_enabled:
            write_split_outputs(merged, resolve_merged_split_output_dir(output_target.path))
        return

    for result in processed_results:
        if args.output_format == "default":
            main_output_path = resolve_raw_default_output_path(
                result.source_path,
                output_target,
                source_label=output_labels[result.source_path],
            )
            write_csv(result.dataframe, main_output_path)
            LOGGER.info("Wrote %s", main_output_path)
            if split_enabled:
                split_target = split_targets[result.source_path]
                write_split_outputs(
                    result.dataframe,
                    split_target.directory,
                    filename_prefix=split_target.filename_prefix,
                )
        elif args.output_format == "taxburst":
            taxburst_dir = resolve_taxburst_output_dir(
                result.source_path,
                output_target,
                source_label=output_labels[result.source_path],
            )
            taxburst_dir.mkdir(parents=True, exist_ok=True)
            main_output_path = resolve_taxburst_main_output_path(
                result.source_path,
                output_target,
                source_label=output_labels[result.source_path],
            )
            write_csv(result.dataframe, main_output_path)
            LOGGER.info("Wrote %s", main_output_path)
            split_target = split_targets[result.source_path]
            write_split_outputs(
                result.dataframe,
                split_target.directory,
                filename_suffix=".taxburst.split.csv",
            )
            LOGGER.info("Wrote taxburst outputs in %s", taxburst_dir)
        else:
            phyloseq_dir = resolve_phyloseq_output_dir(
                result.source_path,
                output_target,
                source_label=output_labels[result.source_path],
            )
            write_phyloseq_outputs(result.dataframe, phyloseq_dir)
            LOGGER.info("Wrote phyloseq outputs in %s", phyloseq_dir)


def process_raw_metaflow_file(
    input_path: Path,
    *,
    use_average: bool,
    remove_unclassified: bool,
    min_coverage: float | None,
    output_format: str = "default",
    suppress_remove_unclassified_warning: bool = False,
) -> ProcessedResult:
    warning_collector = RawMetaflowWarningCollector(input_path)
    dataframe = read_table(input_path)
    validate_raw_metaflow_minimum_columns(dataframe, input_path)
    output_columns = resolve_raw_output_columns(output_format)
    include_synthetic_unclassified = raw_output_includes_synthetic_unclassified(output_format, remove_unclassified)
    required_columns = resolve_raw_metaflow_required_columns(
        use_average=use_average,
        include_synthetic_unclassified=include_synthetic_unclassified,
        min_coverage=min_coverage,
    )
    if min_coverage is not None:
        validate_required_columns(dataframe, ["relative_abundance", "file_name"], input_path)
    warn_on_missing_optional_columns(
        dataframe,
        [column for column in RAW_METAFLOW_REFERENCED_COLUMNS if column not in required_columns],
        input_path,
        context="raw_metaflow",
        collector=warning_collector,
    )
    dataframe = add_missing_columns(dataframe, RAW_METAFLOW_REFERENCED_COLUMNS)
    dataframe = drop_blank_raw_metaflow_records(dataframe, input_path, collector=warning_collector)
    if dataframe.empty:
        standardized = pd.DataFrame(columns=output_columns)
        standardized["__source_sample_name__"] = pd.Series(dtype="object")
        warning_collector.emit()
        return ProcessedResult(
            dataframe=standardized,
            source_path=input_path,
            wgs_ids_complete=True,
        )
    if not wgs_ids_complete(dataframe):
        warning_collector.add(
            "Missing or empty WGS_ID values were found; sample_name will fall back to the input filename where needed."
        )
    dataframe["__sample_name"] = build_sample_name_series(dataframe, input_path)
    dataframe["__provenance"] = str(input_path)
    dataframe = coerce_numeric_columns(
        dataframe,
        RAW_METAFLOW_NUMERIC_COLUMNS,
        input_path,
        allow_missing=True,
    )

    working = dataframe.copy()
    if min_coverage is not None:
        working = apply_min_coverage_filter(working, min_coverage, input_path, collector=warning_collector)

    if working.empty:
        standardized = pd.DataFrame(columns=output_columns)
        standardized["__source_sample_name__"] = pd.Series(dtype="object")
        warning_collector.emit()
        return ProcessedResult(
            dataframe=standardized,
            source_path=input_path,
            wgs_ids_complete=True,
        )

    post_filter_rows = working.copy()
    formula_column = "average_abund" if use_average else "median_abund"
    formula_required_columns = resolve_raw_metaflow_non_null_columns(
        use_average=use_average,
        include_synthetic_unclassified=include_synthetic_unclassified,
        min_coverage=min_coverage,
    )
    working, incomplete_rows = split_formula_ready_rows(working, formula_required_columns)

    if not incomplete_rows.empty:
        if working.empty:
            warning_collector.add(
                f"{len(incomplete_rows)} rows lack formula inputs. Writing passthrough standardized rows with blank computed fields where necessary."
            )
            passthrough = build_passthrough_precomputed_rows(incomplete_rows, output_format=output_format)
            passthrough["__source_sample_name__"] = input_path.stem
            passthrough = passthrough.reindex(columns=output_columns + ["__source_sample_name__"])
            warning_collector.emit()
            return ProcessedResult(
                dataframe=passthrough,
                source_path=input_path,
                wgs_ids_complete=wgs_ids_complete(post_filter_rows),
            )

        warning_collector.add(
            f"{len(incomplete_rows)} rows lack formula inputs and will be skipped from formula-based processing and output."
        )

    require_non_null_values(working, formula_required_columns, input_path)
    ensure_nonzero(working, "scaled", input_path)
    ensure_nonzero(working, "total_weighted_hashes", input_path)
    consistency_group_columns = resolve_raw_metaflow_group_columns(
        working,
        input_path,
        collector=warning_collector,
    )
    ensure_consistent_group_metadata(
        working,
        input_path,
        consistency_group_columns,
        include_synthetic_unclassified=include_synthetic_unclassified,
        collector=warning_collector,
    )
    synthetic_group_columns = resolve_raw_metaflow_group_columns(
        working,
        input_path,
        include_file_name=min_coverage is None,
        warn_on_reduced=False,
        collector=warning_collector,
    )

    working["__computed_relative_abundance"] = (
        (working["unique_intersect_bp"] / working["scaled"]) * working[formula_column]
    ) / working["total_weighted_hashes"]
    working["__n_unique_kmers"] = (working["unique_intersect_bp"] / working["scaled"]) * working[formula_column]

    abundance_for_output = working["__computed_relative_abundance"].copy()
    if remove_unclassified and output_format != "taxburst":
        if not suppress_remove_unclassified_warning:
            warning_collector.add(
                "f_weighted_at_rank has been renormalized because --remove_unclassified was used and must be interpreted carefully."
            )
        abundance_for_output = renormalize_within_group(
            abundance_for_output,
            working["__sample_name"],
            input_path,
            "relative abundance after removing synthetic unclassified rows",
        )

    working["__abundance_for_output"] = abundance_for_output

    standardized = pd.DataFrame(
        {
            "sample_name": working["__sample_name"],
            "lineage": working["lineage"],
            "query_md5": working["query_md5"],
            "query_filename": working["query_filename"],
            "ksize": working["ksize"],
            "scaled": working["scaled"],
            "match_name": working["organism_name"],
            "unweighted_fraction": working["f_unique_to_query"],
            "f_weighted_at_rank": working["__abundance_for_output"],
            "bp_match_at_rank": working["unique_intersect_bp"],
            "total_weighted_hashes": working["total_weighted_hashes"],
            "n_unique_kmers": working["__n_unique_kmers"],
            "p_vals": working["p_vals"],
            "actual_confidence_with_coverage": working["actual_confidence_with_coverage"],
            "alt_confidence_mut_rate_with_coverage": working["alt_confidence_mut_rate_with_coverage"],
            "min_coverage": working["file_name"],
            "source_file": working["__provenance"],
        }
    )
    if output_format == "taxburst":
        standardized["n_unique_weighted_found"] = working["n_unique_weighted_found"]
        standardized["sum_weighted_found"] = working["sum_weighted_found"]
        standardized["median_abund"] = working["median_abund"]
    standardized["__source_sample_name__"] = input_path.stem

    if include_synthetic_unclassified:
        synthetic = build_synthetic_unclassified_rows(
            working,
            input_path=input_path,
            provenance_value=str(input_path),
            blank_min_coverage=min_coverage is not None,
            group_columns=synthetic_group_columns,
            collector=warning_collector,
        )
        synthetic["__source_sample_name__"] = input_path.stem
        standardized = pd.concat([standardized, synthetic], ignore_index=True, sort=False)

    standardized = standardized.reindex(columns=output_columns + ["__source_sample_name__"])
    warning_collector.emit()
    return ProcessedResult(
        dataframe=standardized,
        source_path=input_path,
        wgs_ids_complete=wgs_ids_complete(post_filter_rows),
    )


def apply_min_coverage_filter(
    dataframe: pd.DataFrame,
    min_coverage: float,
    input_path: Path,
    *,
    collector: RawMetaflowWarningCollector | None = None,
) -> pd.DataFrame:
    filtered = dataframe.copy()
    filtered = filtered[filtered["relative_abundance"].notna()].copy()
    filtered = filtered[filtered["file_name"].notna()].copy()

    file_name_text = filtered["file_name"].astype(str)
    token = f"min_coverage{format_value(min_coverage)}"
    token_mask = file_name_text.str.contains(re.escape(token), regex=True, na=False)

    extracted = pd.to_numeric(file_name_text.str.extract(MIN_COVERAGE_PATTERN)[0], errors="coerce")
    numeric_mask = pd.Series(np.isclose(extracted, min_coverage, rtol=0, atol=1e-12), index=filtered.index)

    filtered = filtered[token_mask | numeric_mask].copy()

    if filtered.empty:
        emit_warning(
            f"No rows remained after applying --min_coverage {min_coverage}.",
            collector=collector,
        )
    return filtered


def ensure_nonzero(dataframe: pd.DataFrame, column: str, input_path: Path) -> None:
    invalid = dataframe[column].eq(0)
    if invalid.any():
        raise ProcessingError(
            f"{input_path} contains zero values in '{column}' for {invalid.sum()} rows, which would cause division by zero."
        )


def ensure_consistent_group_metadata(
    dataframe: pd.DataFrame,
    input_path: Path,
    group_columns: Sequence[str],
    *,
    include_synthetic_unclassified: bool = True,
    collector: RawMetaflowWarningCollector | None = None,
) -> None:
    if not group_columns or not include_synthetic_unclassified:
        return

    for column in RAW_METAFLOW_STRICT_GROUP_METADATA_COLUMNS:
        if column not in dataframe.columns or not dataframe[column].notna().any():
            continue
        conflicts = dataframe.groupby(list(group_columns), dropna=False)[column].apply(
            lambda series: series.dropna().nunique() > 1
        )
        if conflicts.any():
            conflict_keys = list(conflicts[conflicts].index[:3])
            raise ProcessingError(
                f"{input_path} contains inconsistent '{column}' values within raw_metaflow grouping units. "
                f"Example conflicting groups: {conflict_keys}"
            )

    for column in RAW_METAFLOW_SOFT_GROUP_METADATA_COLUMNS:
        if column not in dataframe.columns or not dataframe[column].notna().any():
            continue
        conflicts = dataframe.groupby(list(group_columns), dropna=False)[column].apply(
            lambda series: series.dropna().nunique() > 1
        )
        if conflicts.any():
            conflict_keys = list(conflicts[conflicts].index[:3])
            emit_warning(
                f"Inconsistent '{column}' values were found within raw_metaflow grouping units. "
                f"Synthetic unclassified rows leave this metadata blank. Example conflicting groups: {conflict_keys}",
                collector=collector,
            )


def renormalize_within_group(
    values: pd.Series,
    groups: pd.Series,
    input_path: Path,
    context: str,
) -> pd.Series:
    group_sums = values.groupby(groups, dropna=False).transform("sum")
    zero_groups = group_sums.isna() | np.isclose(group_sums, 0.0)
    if zero_groups.any():
        bad_groups = pd.Series(groups[zero_groups]).drop_duplicates().tolist()[:5]
        raise ProcessingError(
            f"{input_path} has {context} group sums equal to zero or missing; cannot renormalize. "
            f"Example groups: {bad_groups}"
        )
    return values / group_sums


def build_synthetic_unclassified_rows(
    working: pd.DataFrame,
    *,
    input_path: Path,
    provenance_value: str,
    blank_min_coverage: bool,
    group_columns: Sequence[str],
    collector: RawMetaflowWarningCollector | None = None,
) -> pd.DataFrame:
    if working.empty:
        return pd.DataFrame(columns=RAW_METAFLOW_STANDARD_COLUMNS)

    grouped = working.groupby(list(group_columns), dropna=False, sort=False)
    summaries = grouped.agg(
        query_bp=("query_bp", first_non_null),
        total_weighted_hashes=("total_weighted_hashes", first_non_null),
        p_vals=("p_vals", first_non_null),
        actual_confidence_with_coverage=("actual_confidence_with_coverage", first_non_null),
        alt_confidence_mut_rate_with_coverage=("alt_confidence_mut_rate_with_coverage", first_non_null),
        sum_unweighted_fraction=("f_unique_to_query", "sum"),
        sum_f_weighted_at_rank=("__abundance_for_output", "sum"),
        sum_bp_match_at_rank=("unique_intersect_bp", "sum"),
    ).reset_index()
    summaries = add_missing_columns(
        summaries,
        [
            "__sample_name",
            "file_name",
            "query_md5",
            "query_filename",
            "ksize",
            "scaled",
            "actual_confidence_with_coverage",
            "alt_confidence_mut_rate_with_coverage",
            "p_vals",
            "total_weighted_hashes",
            "query_bp",
        ],
    )

    summaries["sample_name"] = summaries["__sample_name"]
    summaries["lineage"] = "unclassified"
    summaries["match_name"] = "unclassified"
    summaries["unweighted_fraction"] = 1.0 - summaries["sum_unweighted_fraction"]
    summaries["f_weighted_at_rank"] = 1.0 - summaries["sum_f_weighted_at_rank"]
    summaries["bp_match_at_rank"] = summaries["query_bp"] - summaries["sum_bp_match_at_rank"]
    summaries["n_unique_kmers"] = summaries["f_weighted_at_rank"] * summaries["total_weighted_hashes"]
    summaries["p_vals"] = ""
    summaries["actual_confidence_with_coverage"] = ""
    summaries["alt_confidence_mut_rate_with_coverage"] = ""
    summaries["min_coverage"] = ""
    summaries["source_file"] = provenance_value

    summaries = clip_tiny_negative_values(
        summaries,
        input_path=input_path,
        columns=["unweighted_fraction", "f_weighted_at_rank", "bp_match_at_rank", "n_unique_kmers"],
        collector=collector,
    )
    warn_on_suspicious_unclassified_values(
        summaries,
        input_path=input_path,
        columns=["unweighted_fraction", "f_weighted_at_rank", "bp_match_at_rank", "n_unique_kmers"],
        collector=collector,
    )

    synthetic = summaries[
        [
            "sample_name",
            "lineage",
            "query_md5",
            "query_filename",
            "ksize",
            "scaled",
            "match_name",
            "unweighted_fraction",
            "f_weighted_at_rank",
            "bp_match_at_rank",
            "total_weighted_hashes",
            "n_unique_kmers",
            "p_vals",
            "actual_confidence_with_coverage",
            "alt_confidence_mut_rate_with_coverage",
            "min_coverage",
            "source_file",
        ]
    ].copy()
    return synthetic


def first_non_null(series: pd.Series):
    non_null = series.dropna()
    if non_null.empty:
        return pd.NA
    return non_null.iloc[0]


def clip_tiny_negative_values(
    dataframe: pd.DataFrame,
    *,
    input_path: Path,
    columns: Sequence[str],
    tolerance: float = 1e-9,
    collector: RawMetaflowWarningCollector | None = None,
) -> pd.DataFrame:
    for column in columns:
        if column not in dataframe.columns:
            continue
        mask = dataframe[column].lt(0) & dataframe[column].ge(-tolerance)
        if mask.any():
            emit_warning(
                f"Tiny negative floating-point artifacts were clipped to zero in synthetic unclassified column '{column}' for {int(mask.sum())} rows.",
                collector=collector,
            )
            dataframe.loc[mask, column] = 0.0
    return dataframe


def warn_on_suspicious_unclassified_values(
    dataframe: pd.DataFrame,
    *,
    input_path: Path,
    columns: Sequence[str],
    collector: RawMetaflowWarningCollector | None = None,
) -> None:
    suspicious = pd.Series(False, index=dataframe.index)
    for column in columns:
        if column in dataframe.columns:
            suspicious = suspicious | dataframe[column].lt(0)

    if suspicious.any():
        example_rows = dataframe.loc[suspicious, ["sample_name", "query_md5", "match_name"]].head(5).to_dict("records")
        emit_warning(
            f"Synthetic unclassified values are negative or otherwise suspicious in {int(suspicious.sum())} rows. "
            f"Example rows: {example_rows}",
            collector=collector,
        )


def wgs_ids_complete(dataframe: pd.DataFrame) -> bool:
    if "WGS_ID" not in dataframe.columns:
        return False
    values = dataframe["WGS_ID"].astype("string")
    return values.notna().all() and values.str.strip().ne("").all()


def merge_processed_raw_results(results: Sequence[ProcessedResult]) -> pd.DataFrame:
    ensure_compatible_frames(
        [result.dataframe for result in results],
        context="Processed raw_metaflow outputs are not compatible for merging.",
    )

    merged = pd.concat([result.dataframe for result in results], ignore_index=True, sort=False)
    use_wgs_id = all(result.wgs_ids_complete for result in results)
    if not use_wgs_id:
        merged["sample_name"] = merged["__source_sample_name__"]

    output_columns = [column for column in merged.columns if column != "__source_sample_name__"]
    merged = merged.reindex(columns=output_columns + ["__source_sample_name__"])
    merged = merged.drop(columns="__source_sample_name__")
    return merged


def ensure_compatible_frames(frames: Sequence[pd.DataFrame], *, context: str) -> None:
    if not frames:
        return
    expected = list(frames[0].columns)
    for index, frame in enumerate(frames[1:], start=2):
        if list(frame.columns) != expected:
            raise ProcessingError(f"{context} Frame {index} has a different column order/schema.")


def build_raw_output_labels(source_paths: Sequence[Path], output_target: OutputTarget) -> dict[Path, str]:
    labels = {path: path.name for path in source_paths}
    if output_target.mode != "directory":
        return labels

    grouped_by_name: dict[str, list[Path]] = {}
    for source_path in source_paths:
        grouped_by_name.setdefault(source_path.name, []).append(source_path)

    for name, paths in grouped_by_name.items():
        if len(paths) == 1:
            continue
        for path in paths:
            labels[path] = f"{path.stem}__{short_path_hash(path)}{path.suffix}"
    return labels


def resolve_raw_default_output_path(
    source_path: Path,
    output_target: OutputTarget,
    *,
    source_label: str | None = None,
) -> Path:
    label = source_label or source_path.name
    filename = f"{label}.processed_metaflow.csv"
    if output_target.mode == "auto":
        return source_path.parent / filename
    if output_target.mode == "directory":
        assert output_target.path is not None
        return output_target.path / filename
    assert output_target.path is not None
    return output_target.path


def resolve_phyloseq_output_dir(
    source_path: Path,
    output_target: OutputTarget,
    *,
    source_label: str | None = None,
) -> Path:
    label = source_label or source_path.name
    directory_name = f"{label}.phyloseq"
    if output_target.mode == "auto":
        return source_path.parent / directory_name
    assert output_target.path is not None
    return output_target.path / directory_name


def resolve_taxburst_output_dir(
    source_path: Path,
    output_target: OutputTarget,
    *,
    source_label: str | None = None,
) -> Path:
    label = source_label or source_path.name
    directory_name = f"{label}.taxburst"
    if output_target.mode == "auto":
        return source_path.parent / directory_name
    assert output_target.path is not None
    return output_target.path / directory_name


def resolve_taxburst_main_output_path(
    source_path: Path,
    output_target: OutputTarget,
    *,
    source_label: str | None = None,
) -> Path:
    label = source_label or source_path.name
    return resolve_taxburst_output_dir(source_path, output_target, source_label=source_label) / f"{label}.taxburst.csv"


def build_split_targets(
    processed_results: Sequence[ProcessedResult],
    output_target: OutputTarget,
    output_format: str,
    output_labels: dict[Path, str],
) -> dict[Path, SplitTarget]:
    if output_format not in {"default", "taxburst"}:
        return {}

    if output_format == "taxburst":
        return {
            result.source_path: SplitTarget(
                directory=resolve_taxburst_output_dir(
                    result.source_path,
                    output_target,
                    source_label=output_labels[result.source_path],
                )
            )
            for result in processed_results
        }

    split_directories: dict[Path, Path] = {
        result.source_path: resolve_split_output_dir(result.source_path, output_target)
        for result in processed_results
    }
    grouped_sources: dict[Path, list[Path]] = {}
    for source_path, directory in split_directories.items():
        grouped_sources.setdefault(directory, []).append(source_path)

    split_targets: dict[Path, SplitTarget] = {}
    for directory, source_paths in grouped_sources.items():
        prefixes = build_split_prefixes(source_paths) if len(source_paths) > 1 else {source_paths[0]: ""}
        for source_path in source_paths:
            split_targets[source_path] = SplitTarget(
                directory=directory,
                filename_prefix=prefixes[source_path],
            )
    return split_targets


def resolve_split_output_dir(source_path: Path, output_target: OutputTarget) -> Path:
    if output_target.mode == "auto":
        return source_path.parent / "splits"
    assert output_target.path is not None
    if output_target.mode == "directory":
        return output_target.path / "splits"
    return output_target.path.parent / "splits"


def build_split_prefixes(source_paths: Sequence[Path]) -> dict[Path, str]:
    grouped_by_label: dict[str, list[Path]] = {}
    for source_path in source_paths:
        label = sanitize_filename(source_path.stem) or "input"
        grouped_by_label.setdefault(label, []).append(source_path)

    prefixes: dict[Path, str] = {}
    for label, paths in grouped_by_label.items():
        if len(paths) == 1:
            prefixes[paths[0]] = f"{label}__"
            continue
        for path in paths:
            prefixes[path] = f"{label}__{short_path_hash(path)}__"
    return prefixes


def resolve_merged_split_output_dir(output_file: Path) -> Path:
    return output_file.parent / "splits"


def write_split_outputs(
    dataframe: pd.DataFrame,
    output_directory: Path,
    *,
    filename_prefix: str = "",
    filename_suffix: str = ".processed_metaflow.split.csv",
) -> None:
    if "__source_sample_name__" in dataframe.columns:
        frame_to_split = dataframe.drop(columns="__source_sample_name__")
    else:
        frame_to_split = dataframe

    output_paths: list[Path] = []
    payloads: list[tuple[Path, pd.DataFrame]] = []
    sample_series = frame_to_split["sample_name"].astype("object")
    unique_sample_names = list(pd.unique(sample_series))

    for sample_name in unique_sample_names:
        if pd.isna(sample_name):
            group_frame = frame_to_split[sample_series.isna()].copy()
        else:
            group_frame = frame_to_split[sample_series == sample_name].copy()
        label = sanitize_filename(format_value(sample_name))
        if not label:
            label = "unnamed_sample"
        output_path = output_directory / f"{filename_prefix}{label}{filename_suffix}"
        output_paths.append(output_path)
        payloads.append((output_path, group_frame.copy()))

    ensure_unique_output_paths(output_paths)
    output_directory.mkdir(parents=True, exist_ok=True)

    for output_path, split_frame in payloads:
        write_csv(split_frame, output_path)
        LOGGER.info("Wrote %s", output_path)


def write_phyloseq_outputs(dataframe: pd.DataFrame, output_dir: Path) -> None:
    if "__source_sample_name__" in dataframe.columns:
        dataframe = dataframe.drop(columns="__source_sample_name__")

    validate_phyloseq_input(dataframe)
    output_dir.mkdir(parents=True, exist_ok=True)

    tax_table = build_tax_table(dataframe)
    otu_table = build_otu_table(dataframe)
    sample_data = build_sample_data(dataframe)

    write_csv(tax_table, output_dir / "tax_table.csv")
    write_csv(otu_table, output_dir / "otu_table.csv")
    write_csv(sample_data, output_dir / "sample_data.csv")


def validate_phyloseq_input(dataframe: pd.DataFrame) -> None:
    missing_lineage = dataframe["lineage"].isna() | dataframe["lineage"].astype(str).str.strip().eq("")
    if missing_lineage.any():
        raise ProcessingError(
            "Cannot create phyloseq output because at least one processed record has missing lineage. "
            "Inspect the input data before retrying."
        )


def build_tax_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    distinct = dataframe[["match_name", "lineage"]].drop_duplicates().copy()
    conflicts = distinct.groupby("match_name", dropna=False)["lineage"].nunique(dropna=False).gt(1)
    if conflicts.any():
        taxa = list(conflicts[conflicts].index[:5])
        raise ProcessingError(
            f"Cannot build phyloseq tax_table because match_name values map to multiple lineages: {taxa}"
        )

    split_lineage = distinct["lineage"].astype(str).str.split(";", expand=True)
    if split_lineage.shape[1] > 7:
        raise ProcessingError("Lineage values contain more than seven ranks and cannot be mapped to phyloseq columns.")

    split_lineage = split_lineage.reindex(columns=range(7), fill_value="")
    split_lineage.columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    tax_table = pd.concat([distinct[["match_name"]].reset_index(drop=True), split_lineage.reset_index(drop=True)], axis=1)
    return tax_table


def build_otu_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    duplicates = dataframe.duplicated(subset=["match_name", "sample_name"], keep=False)
    if duplicates.any():
        LOGGER.warning(
            "Multiple records share the same match_name/sample_name combination; phyloseq otu_table will sum n_unique_kmers within each combination."
        )

    otu = (
        dataframe.pivot_table(
            index="match_name",
            columns="sample_name",
            values="n_unique_kmers",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    return otu


def build_sample_data(dataframe: pd.DataFrame) -> pd.DataFrame:
    columns = ["sample_name"]
    if "source_file" in dataframe.columns:
        columns.append("source_file")
    return dataframe[columns].drop_duplicates().reset_index(drop=True)


def ensure_unique_output_paths(paths: Iterable[Path]) -> None:
    seen: set[Path] = set()
    duplicates: list[Path] = []
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            duplicates.append(resolved)
        seen.add(resolved)
    if duplicates:
        raise ProcessingError(
            f"Output path collision detected. Repeated destinations include: {', '.join(str(path) for path in duplicates[:5])}"
        )


def write_csv(dataframe: pd.DataFrame, path: Path) -> None:
    output = dataframe.copy()
    if "__source_sample_name__" in output.columns:
        output = output.drop(columns="__source_sample_name__")
    output = normalize_csv_output(output)
    output.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")


def normalize_csv_output(dataframe: pd.DataFrame) -> pd.DataFrame:
    normalized = dataframe.copy()
    for column in normalized.columns:
        normalized[column] = normalized[column].map(format_value)
    return normalized


def sanitize_filename(value: str) -> str:
    sanitized = INVALID_FILENAME_CHARS.sub("_", value).strip().rstrip(".")
    sanitized = re.sub(r"\s+", "_", sanitized)
    return sanitized or "unnamed"


def short_path_hash(path: Path) -> str:
    return hashlib.sha1(str(path.resolve()).encode("utf-8")).hexdigest()[:8]


def format_value(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if float(value).is_integer():
            return str(int(value))
        return format(float(value), ".15g")
    text = str(value).strip()
    if re.fullmatch(r"-?\d+\.0+", text):
        return text.split(".", maxsplit=1)[0]
    return text


if __name__ == "__main__":
    raise SystemExit(main())
