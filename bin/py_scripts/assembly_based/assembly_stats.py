#!/usr/bin/env python

# Calculates assembly statistics for each assembled genome in a directory.

__version__ = '0.1.1'
__date__ = '27-10-2025'
__author__ = 'D.J.BERGER, Quan TNM'

###### Imports

import argparse
import csv
import gzip
import sys
from itertools import chain
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

###### Functions

# Scaffold stats
def calculate_scaffold_stats(records: list) -> tuple:
    """
    Calculate scaffold statistics from parsed records.
    """
    count = len(records)
    total_asm = sum(len(record.seq) for record in records)
    lengths = sorted([len(record.seq) for record in records], reverse=True)
    gc_count = sum(
        record.seq.upper().count("G") + record.seq.upper().count("C")
        for record in records
    )
    all_count = sum(
        record.seq.upper().count("G") + record.seq.upper().count("C") +
        record.seq.upper().count("T") + record.seq.upper().count("A")
        for record in records
    )
    gc_cont = round((gc_count / all_count) * 100, 2) if all_count > 0 else 0.0
    # Avoid division by zero
    s_n50 = (
        next(
            (length for length in lengths
             if sum(lengths[:lengths.index(length) + 1]) >= total_asm * 0.5),
            None
        ) if lengths else None
    )
    s_n90 = (
        next(
            (length for length in lengths
             if sum(lengths[:lengths.index(length) + 1]) >= total_asm * 0.9),
            None
        ) if lengths else None
    )
    gap_sum = sum(record.seq.count("N") for record in records)
    gap_perc = round((gap_sum / total_asm) * 100, 2) if total_asm > 0 else 0.0
    return count, s_n50, gap_sum, total_asm, s_n90, gc_cont, gap_perc

# Contig stats
def calculate_contig_stats(records: list, gap_size: int) -> tuple:
    """
    Calculate contig statistics by splitting scaffolds on gaps.
    Optimized to avoid creating new records; compute lengths directly.
    """
    contig_lengths = []
    gap_count = 0
    for record in records:
        seq = str(record.seq)
        contigs = seq.split("N" * gap_size)
        for contig in contigs:
            if len(contig) > 0 and "N" not in contig:
                contig_lengths.append(len(contig))
        # Count gaps as splits beyond original scaffolds
        gap_count += len(contigs) - 1
    c_count = len(contig_lengths)
    contig_lengths.sort(reverse=True)
    total_len = sum(contig_lengths)
    c_n50 = (
        next(
            (length for length in contig_lengths
             if sum(contig_lengths[:contig_lengths.index(length) + 1]) >= total_len * 0.5),
            None
        ) if contig_lengths else None
    )
    c_n90 = (
        next(
            (length for length in contig_lengths
             if sum(contig_lengths[:contig_lengths.index(length) + 1]) >= total_len * 0.9),
            None
        ) if contig_lengths else None
    )
    return c_count, c_n50, gap_count, c_n90

# Write results to csv file
def write_stats_to_file(
    bin_name: str,
    count: int,
    s_n50,
    c_count: int,
    c_n50,
    gap_count: int,
    gap_sum: int,
    total_asm: int,
    s_n90,
    c_n90,
    gc_cont: float,
    gap_perc: float,
    csvfile
):
    writer = csv.writer(csvfile)
    writer.writerow([
        bin_name, total_asm, count, s_n50, s_n90, c_count, c_n50, c_n90,
        gc_cont, gap_count, gap_sum, gap_perc
    ])

# Parse arguments from the command line.
def parse_args():
    description = (
        f'Calculates assembly statistics for each assembled genome in a directory. '
        f'Version: {__version__}, Date: {__date__}, Author: {__author__}'
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '--fasta_dir', required=True, help='Directory containing FASTA files'
    )
    parser.add_argument(
        '--gap', type=int, default=2,
        help="Minimum gap length to be considered a scaffold [2]"
    )
    parser.add_argument(
        '--output', default="sample", help="Output file prefix ['sample']"
    )
    parser.add_argument(
        "--version", action="version", version=f'Version: {__version__}'
    )
    return parser.parse_args()

###### Main
def main():
    args = parse_args()
    fasta_dir = Path(args.fasta_dir)
    if not fasta_dir.is_dir():
        sys.exit(f"Error: {fasta_dir} is not a valid directory.")

    with open(f"{args.output}.assembly_stats.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "bin_name", "assembly_length_bp", "scaffold_count", "scaffold_N50_bp",
            "scaffold_N90_bp", "contig_count", "contig_N50_bp", "contig_N90_bp",
            "GC_perc", "gaps_count", "gaps_sum_bp", "gaps_perc"
        ])

        # Find and process FASTA files (supporting .fasta, .fa, .fna, and their .gz variants)
        fasta_patterns = [
            "*.fasta", "*.fa", "*.fna", "*.fasta.gz", "*.fa.gz", "*.fna.gz"
        ]
        for pattern in fasta_patterns:
            for file_path in fasta_dir.glob(pattern):
                try:
                    # Extract bin_name as the filename without extensions
                    bin_name = Path(file_path.stem).stem
                    # Parse records, handling gzipped files
                    if str(file_path).endswith('.gz'):
                        with gzip.open(file_path, 'rt') as f:
                            records = list(SeqIO.parse(f, "fasta"))
                    else:
                        records = list(SeqIO.parse(file_path, "fasta"))
                    if not records:
                        print(f"Warning: No records in {file_path.name}. Skipping.")
                        continue
                    count, s_n50, gap_sum, total_asm, s_n90, gc_cont, gap_perc = (
                        calculate_scaffold_stats(records)
                    )
                    c_count, c_n50, gap_count, c_n90 = (
                        calculate_contig_stats(records, args.gap)
                    )
                    write_stats_to_file(
                        bin_name, count, s_n50, c_count, c_n50, gap_count,
                        gap_sum, total_asm, s_n90, c_n90, gc_cont, gap_perc, csvfile
                    )
                except Exception as e:
                    print(f"Error processing {file_path.name}: {e}. Skipping.")

if __name__ == "__main__":
    main()