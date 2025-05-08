# profiler_sourmash

A read-based shotgun metagenomics processing pipeline built using Nextflow DSL2, combining sourmash and YACHT tools for efficient taxonomic profiling.

## Overview

This pipeline performs taxonomic profiling of shotgun metagenomics data through the following steps:
1. Merging paired-end sequences
2. Generating sketches using sourmash
3. Performing fast multi-gather analysis with sourmash
4. Running YACHT analysis
5. Processing and combining results

## Quick Start

```bash
# Using CSV input
nextflow run main.nf \
    --input samples.csv \
    --input_format csv \
    -profile conda

# Using directory input
nextflow run main.nf \
    --trimmed_fastq /path/to/fastq/directory \
    --input_format directory \
    -profile conda
```

## Requirements
* Nextflow (>=23.04.0)
* One of the following container engines:
    * Conda (with mamba)
    * Singularity
    * Docker
    * Apptainer
## Input Requirements
* Sample CSV format:
```
sample_id,run_id,group,short_reads_1,short_reads_2,long_reads
sample1,1,0,data/sample1_1.fq.gz,data/sample1_2.fq.gz,data/sample1.fastq.gz
sample1,2,0,data/sample1_1.fastq.gz,data/sample1_2.fastq.gz,data/sample1.fastq.gz
sample2,0,0,data/sample2_1.fastq.gz,data/sample2_2.fastq.gz,data/sample2.fastq.gz
sample3,1,1,data/sample3_1.fastq.gz,data/sample3_2.fastq.gz,
```
* Diectory containing fastq files.

### Supported File Formats
* Short reads: `.fastq.gz`, `.fastq`
* Long reads: `.fasta`, `.fa`

## Pipeline Parameters

### Required Parameters

| Parameter                  | Description                               | Default Value |
| -------------------------- | ----------------------------------------- | ------------- |
| `input` or `trimmed_fastq` | Path to input CSV file or FastQ directory | `null`        |
| `input_format`             | Format of input (`csv` or `directory`)    | `csv`         |
| `sourmash_database`        | Path to sourmash RocksDB database         | `null`        |
| `yacht_database`           | Path to YACHT database JSON               | `null`        |

### Directory Configuration

| Parameter          | Description                         | Default Value                           |
| ------------------ | ----------------------------------- | --------------------------------------- |
| `base_dir`         | Base directory for pipeline outputs | `null`                                  |
| `workDir`          | Working directory for pipeline      | `null`                                  |
| `outdir`           | Output directory                    | `${params.base_dir}/results`            |
| `results_dir`      | Directory for results               | `${params.base_dir}/results`            |
| `intermediate_dir` | Directory for intermediate files    | `${params.base_dir}/intermediate_files` |
| `sketch_dir`       | Directory for sketch files          | `${params.base_dir}/sketches`           |
| `merged_seq_dir`   | Directory for merged sequences      | `${params.base_dir}/merged_seq`         |

### Analysis Parameters

| Parameter | Description             | Default Value |
| --------- | ----------------------- | ------------- |
| `ksize`   | K-mer size for sourmash | `51`          |

### Process Control

| Parameter                 | Description                      | Default Value |
| ------------------------- | -------------------------------- | ------------- |
| `enable_copymergedseqs`   | Copy merged sequence files       | `true`        |
| `enable_copysketch`       | Copy sketch files                | `true`        |
| `enable_copyintermediate` | Copy intermediate files          | `true`        |
| `publish_dir_mode`        | Mode for publishing output files | `copy`        |

### Cleanup Options

| Parameter          | Description                          | Default Value                                                               |
| ------------------ | ------------------------------------ | --------------------------------------------------------------------------- |
| `cleanup`          | Enable cleanup of temporary files    | `false`                                                                     |
| `cleanup_patterns` | Patterns to clean up                 | `['work', 'intermediate_files', '.nextflow*', '.command.*', 'process_tmp']` |
| `cleanup_dry_run`  | Preview cleanup without deleting     | `false`                                                                     |
| `keep_days`        | Keep files newer than specified days | `7`                                                                         |

### Notification

| Parameter | Description             | Default Value |
| --------- | ----------------------- | ------------- |
| `email`   | Email for notifications | `null`        |

## Pipeline Profiles

- `conda`: Uses conda/mamba
- `singularity`: Uses Singularity containers
- `mamba`: Uses mamba for faster conda setup
- `apptainer`: Uses Apptainer containers
- `standard`: Default profile (Singularity-enabled)

## Output Structure
```
results/
├── merged_seq/          # Merged paired-end sequences
├── sketches/           # Sourmash sketch files
├── intermediate_files/ # Processing intermediates
│   └── yacht/         # YACHT analysis results
└── processed_results/  # Final analysis results
```

## Resource Requirements

Default computational resources per process:

- Low: 2 CPUs, 4GB memory
- Medium: 16 CPUs, 32GB memory
- High: 20 CPUs, 60GB memory

## Error Handling

The pipeline implements automatic error recovery:

- Maximum 3 retries per task
- Automatic resource scaling on retry
- Configurable error strategy

## Citation

If you use this pipeline, please cite:

## Author
* Ton Ngoc Minh Quan
* Mai Thu Si Nguyen