# `sourmash` Profiler

A read-based shotgun metagenomics processing pipeline built using Nextflow DSL2, combining sourmash and YACHT tools for efficient taxonomic profiling.

## Overview

This pipeline performs taxonomic profiling of shotgun metagenomics data through the following stages:
1.  **Preprocessing:**
    *   Raw read quality control (FastQC).
    *   Adapter and quality trimming (fastp).
    *   Host genome removal (hostile).
    *   Quality control on cleaned reads.
2.  **Sequence Merging:** Merges paired-end reads.
3.  **Sketch Generation:** Creates sequence sketches using `sourmash`.
4.  **Gather Analysis:** Performs fast multi-gather analysis with `sourmash`.
5.  **YACHT Analysis:** Runs taxonomic classification with `YACHT`.
6.  **Result Processing:** Aggregates and finalizes analysis results.

## Quick Start

```bash
# With csv file
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --hostile_reference /path/to/hostile/reference \
    --hostile_index "your-hostile-index-name" \
	--hostile_ref_name "your-hostile-ref-name" \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
	--ksize 31

# With directory of FASTQ files
nextflow run main.nf \
    --input /path/to/your/fastq_directory \
    --input_format directory \
    -profile conda \
    --outdir /path/to/output/directory \
    --hostile_reference /path/to/hostile/reference \
    --hostile_index "your-hostile-index-name" \
    --hostile_ref_name "your-hostile-ref-name" \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --ksize 31
```

## Requirements
*   Nextflow (>=23.04.0)
*   One of the following container engines:
    *   Conda (with mamba)
    *   Singularity
    *   Docker
    *   Apptainer

## Input Specifications
*   **Sample CSV Format:**
    ```csv
    sample_id,run_id,group,short_reads_1,short_reads_2,long_reads
    sample1,1,0,data/sample1_1.fq.gz,data/sample1_2.fq.gz,data/sample1.fastq.gz
    sample2,0,0,data/sample2_1.fastq.gz,data/sample2_2.fastq.gz,
    ```
*   **Directory:** A directory containing paired-end FASTQ files.

### Supported File Formats
*   Short reads: `.fastq.gz`, `.fastq`
*   Long reads: `.fasta`, `.fa`

## Pipeline Parameters

### Required Parameters

| Parameter           | Description                               | Default Value |
| ------------------- | ----------------------------------------- | ------------- |
| `input`             | Path to input CSV file or FASTQ directory | `null`        |
| `input_format`      | Input format (`csv` or `directory`)       | `csv`         |
| `sourmash_database` | Path to sourmash SBT database             | `null`        |
| `yacht_database`    | Path to YACHT database JSON file          | `null`        |

### Preprocessing Parameters

| Parameter             | Description                                                        | Default Value     |
| --------------------- | ------------------------------------------------------------------ | ----------------- |
| `hostile_reference`   | Path to a pre-existing hostile reference directory. If `null`, it will be downloaded. | `null`            |
| `hostile_index`       | Name of the hostile index to use or download.                      | `"human-t2t-hla"` |
| `hostile_ref_name`    | Explicit reference name for hostile.                               | `"human-t2t-hla"` |
| `fastp_adapter_fasta` | Path to a FASTA file containing adapter sequences for fastp.       | `null`            |

### Directory Configuration

| Parameter               | Description                               | Default Value                               |
| ----------------------- | ----------------------------------------- | ------------------------------------------- |
| `outdir`                | Main output directory                     | `./results`                                 |
| `merged_seq_dir`        | Directory for merged sequences            | `${params.outdir}/Merged sequences`         |
| `sketches_dir`          | Directory for sourmash sketch files       | `${params.outdir}/Sourmash - YACHT/sketches`|
| `yacht_results_dir`     | Directory for YACHT results               | `${params.outdir}/Sourmash - YACHT`         |
| `fastmultigather_dir`   | Directory for fastmultigather results     | `${params.outdir}/Sourmash - YACHT`         |
| `processed_results_dir` | Directory for final processed results     | `${params.outdir}/Sourmash - YACHT`         |
| `qc_raw_dir`            | Directory for raw read QC reports         | `${params.outdir}/QC/Raw_reads`             |
| `qc_trim_dir`           | Directory for trimmed read QC reports     | `${params.outdir}/QC/Trimming`              |
| `qc_rmhost_dir`         | Directory for host-removed QC reports     | `${params.outdir}/QC/Remove host genome`    |
| `trim_dir`              | Directory for trimmed reads               | `${params.outdir}/Trimming`                 |
| `rmhost_dir`            | Directory for host-removed reads          | `${params.outdir}/Remove host genome`       |

### Analysis Parameters

| Parameter | Description             | Default Value |
| --------- | ----------------------- | ------------- |
| `ksize`   | K-mer size for sourmash | `31`          |

### Process Control

| Parameter               | Description                      | Default Value |
| ----------------------- | -------------------------------- | ------------- |
| `enable_copymergedseqs` | Copy merged sequence files       | `true`        |
| `enable_copysketch`     | Copy sketch files                | `true`        |
| `enable_copyintermediate` | Copy intermediate files          | `true`        |
| `publish_dir_mode`      | Mode for publishing output files | `copy`        |

### Cleanup Options

| Parameter          | Description                          | Default Value                                                               |
| ------------------ | ------------------------------------ | --------------------------------------------------------------------------- |
| `cleanup`          | Enable cleanup of temporary files    | `false`                                                                     |
| `cleanup_patterns` | Glob patterns for files to clean up  | `['work', 'intermediate_files', '.nextflow*', '.command.*', 'process_tmp']` |
| `cleanup_dry_run`  | Preview cleanup without deleting     | `false`                                                                     |
| `keep_days`        | Keep files newer than specified days | `7`                                                                         |

### Notification

| Parameter | Description             | Default Value |
| --------- | ----------------------- | ------------- |
| `email`   | Email for notifications | `null`        |

### Customizing Tool Arguments
For fine-grained control over the arguments passed to specific tools, define them in a Nextflow configuration file. Use the `process` scope, the `withName` selector, and the `ext` map to pass custom arguments.

For example, create a custom configuration file and supply it to the run command:
```bash
nextflow run main.nf \
    --input samples.csv \
    -profile conda \
    -c /path/to/your/custom_params.config
```

## Pipeline Profiles

-   `standard`: Default profile, uses Singularity.
-   `conda`: Uses conda/mamba with the provided environment file.
-   `singularity`: Uses Singularity containers.
-   `mamba`: Uses mamba for faster conda environment setup.
-   `apptainer`: Uses Apptainer containers.
-   `docker`: Uses Docker containers.

## Output Structure
```
<outdir>/
├── Merged sequences/
│   └── [sample_id]/
│       └── *.merged.fastq
├── QC/
│   ├── Raw_reads/
│   │   └── [sample_id]/
│   │       └── *fastqc.{html,zip}
│   ├── Remove host genome/
│   │   └── [sample_id]/
│   │       └── *fastqc.{html,zip}
│   └── Trimming/
│       └── [sample_id]/
│           └── *fastqc.{html,zip}
├── Remove host genome/
│   └── [sample_id]/
│       ├── cleaned_reads/
│       │   └── *.fastq.gz
│       └── report/
│           └── *.json
├── Sourmash - YACHT/
│   ├── manysketch_manifest.csv
│   ├── sketches/
│   │   └── *.sig
│   ├── [sample_id]_fastmultigather_results.csv
│   ├── [sample_id]_yacht_results.xlsx
│   └── Final_processed_results.xlsx
└── Trimming/
    └── [sample_id]/
        ├── *.fastq.gz
        └── report/
            └── *.{html,json}
```

## Resource Requirements

Default compute resources are assigned per process via labels. These values may increase automatically on retry.

-   `process_low`: 4 CPUs, 12GB memory
-   `process_medium`: 12 CPUs, 60GB memory
-   `process_high`: 25 CPUs, 120GB memory

Maximum limits can be set via the `max_cpus`, `max_memory`, and `max_time` parameters.

## Error Handling

The pipeline implements an automatic error recovery strategy:

-   Maximum of 3 retries per task.
-   Automatic resource scaling on retry.
-   Configurable error strategy.

## Citation

If you use this pipeline, please cite:

## Author
*   Ton Ngoc Minh Quan
*   Mai Thu Si Nguyen