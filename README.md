# metaflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with slurm](https://img.shields.io/badge/run%20with-slurm-1AAEE8.svg?labelColor=000000)](https://www.schedmd.com)
---

**metaflow** is a robust, modular pipeline built on [Nextflow](https://www.nextflow.io/) for comprehensive analysis of Illumina **short-read** metagenomic data. Designed with practical experience in mind, this workflow streamlines both rapid, read-level classification and high-resolution genome assembly, making it accessible to users with minimal computational background. It delivers fast, actionable results for professionals in molecular epidemiology, enabling timely insights for outbreak investigation and surveillance.

This pipeline was developed based on the concepts and design principles of [`nf-core/mag`](https://nf-co.re/mag/latest/) and [`SOMA`](https://github.com/ukhsa-collaboration/SOMA). Many modifications have been made to the original pipelines, reflecting both practical experience and methodological preferences. Most importantly, **metaflow** is optimized to run efficiently on resource-limited devices, ensuring broad accessibility in both laboratory and field settings.


<p align="center">
  <img src="https://github.com/tnmquann/profiler_sourmash/blob/main/docs/images/metaflow_currentbeta_cleaned_metromap_light.png" width="70%" alt="metaflow workflow diagram"/>
</p>

## Overview

This pipeline performs taxonomic profiling of shotgun metagenomics data through the following stages:
- **Read-based Analysis**  
  - Taxonomic assignment for each sequencing read is performed using an optimized combination of [`sourmash`](https://github.com/sourmash-bio/sourmash) and [`YACHT`](https://github.com/KoslickiLab/YACHT). This approach ensures fast, scalable processing, even for large datasets. Reference databases can be customized to fit your specific research needs.
  - Detection of antimicrobial resistance genes (ARGs) is integrated via [`rgi_bwt`](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst), providing actionable insights into resistance profiles.

- **Metagenome Assembly & Binning**  
  - The pipeline assembles contigs and performs binning to recover Metagenome Assembled Genomes (MAGs), enabling *in silico* characterization of microbial communities.

For step-by-step installation instructions, parameter explanations, and advanced usage tutorials, please refer to [wiki page](https://github.com/tnmquann/profiler_sourmash/wiki).

## Pipeline requirements and Installation overview

### Software environment
- Nextflow version 24 or later is required as the workflow engine and depends on Java versions 16–22 for execution stability.
- Conda or Mamba is used for reproducible package and environment management, with Mamba preferred for faster dependency resolution.
- Container runtimes such as Docker or Apptainer are optional but improve portability and reproducibility across heterogeneous systems.

### Hardware and System requirements
- A POSIX-compatible operating system is required, with Linux or macOS preferred; Windows is supported via WSL for compatibility.
- Memory requirements scale with the assembler: a minimum of 32 GB RAM for MEGAHIT and up to 128 GB for metaSPAdes.
- Disk usage is substantial, requiring at least 256 GB to accommodate databases, environments, intermediate files, and per-sample outputs; GPU acceleration is optional for select assembly-based steps.

### Installation workflow
- The pipeline can be obtained via Git cloning, downloading a release archive, or pulling directly with Nextflow for version-controlled deployment.
- Nextflow installation involves validating Java, downloading the executable, optionally adding it to the system PATH, and verifying functionality.
- Conda or Mamba installation finalizes the environment setup, with optional container runtime installation to enhance execution consistency.
> [!NOTE]
> metaflow relies on `sourmash` and `YACHT` for robust taxonomic classification. For detailed setup instructions, please refer to the [Manual Database Setup](https://github.com/tnmquann/profiler_sourmash/wiki/Installation#manual-database-setup-instructions) section.

## Input specifications
*   **Sample CSV format:**
    ```csv
    sample_id,run_id,group,short_reads_1,short_reads_2,long_reads
    sample1,1,0,data/sample1_1.fq.gz,data/sample1_2.fq.gz,data/sample1.fastq.gz
    sample2,0,0,data/sample2_1.fastq.gz,data/sample2_2.fastq.gz,
    ```
*   **Directory:** A directory containing paired-end FASTQ files.

## Quick start

### Read-based analysis
Once your `sourmash` and `YACHT` databases are ready, you can launch the `metaflow` pipeline. Below are examples for different input formats.
```bash
# With csv file
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --enable_readbase

# With directory of FASTQ files
nextflow run main.nf \
    --input /path/to/your/fastq_directory \
    --input_format directory \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --enable_readbase
```

### Assembly-based analysis
For assembly-based workflows, ensure your `sourmash` database is prepared.
```bash
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database
```

## Output structure
### Read-based subworkflow
The following tree illustrates the organization of output files generated by the read-based subworkflow:

```text
outdir/
├── Databases/
│   └── RGI/
├── Preprocess/
│   ├── Merged_sequences/
│   └── QC/
│       ├── Raw_reads/
│       ├── Remove_HostGenome/
│       └── Trimming/
├── Remove_HostGenome/
├── Trimming/
├── Readbased_Analysis/
│   ├── rgi_bwt/
│   └── Sourmash-YACHT/
│       ├── final_results/
│       ├── fastmultigather/
│       ├── reads_collected/
│       ├── sketches/
│       │   ├── single_sketches/
│       │   └── batch.batch.manysketch.zip
│       ├── taxannotate/
│       ├── taxmetagenome/
│       └── yacht_results/
```

Each directory is created automatically by the pipeline and contains outputs relevant to its analysis step. The results of the ARG read-based analysis will be located in the folder `./Readbased_Analysis/rgi_bwt/`, and the results of the read-based taxonomic classification will be located in `./Readbased_Analysis/Sourmash-YACHT/final_results/`.

### Assembly-based subworkflow

The following tree illustrates the organization of output files generated by the assembly-based subworkflow:

```text
outdir/
├── Databases/
│   ├── genomad/
│   ├── bakta/
│   ├── CheckM/
│   ├── CheckM2/
│   ├── GUNC/
│   └── BUSCO/
├── Preprocess/
│   ├── Merged_sequences/
│   └── QC/
│       ├── Raw_reads/
│       ├── Remove_HostGenome/
│       └── Trimming/
├── Remove_HostGenome/
├── Trimming/
├── Assembly/
│   ├── MEGAHIT/
│   │   └── QC/
│   └── SPAdes/
│       └── QC/
├── Binning/
│   ├── Annotation/
│   │   ├── Prokka/
│   │   └── Bakta/
│   ├── Binette/
│   ├── DASTool/
│   ├── COMEBin/
│   ├── MaxBin2/
│   ├── MetaBAT2/
│   ├── SemiBin2/
│   ├── CONCOCT/
│   ├── VAMB/
│   ├── depths/
│   ├── mapping/
│   └── QC/
│       ├── AssemblyStats/
│       ├── CheckM/
│       ├── CheckM2/
│       ├── QUAST/
│       ├── GUNC/
│       └── BUSCO/
├── ContigsAnalysis/
│   ├── annotation/
│   │   ├── Prokka/
│   │   └── pyrodigal/
│   ├── contig_stats/
│   ├── coverage/
│   ├── genomad/
│   ├── skani/
│   └── mapped_reads/
├── Taxonomy/
│   ├── sourmash/
│   └── Tiara/
```

## Error Handling

The pipeline implements an automatic error recovery strategy:
-   Maximum of 3 retries per task.
-   Automatic resource scaling on retry.
-   Configurable error strategy.

## Author
*   Ton Ngoc Minh Quan
*   Mai Thu Si Nguyen

## Citation

If you use this pipeline, please cite:
> Minh-Quan T.N., Nguyen M.T.S. (2025). metaflow: A Nextflow pipeline for comprehensive analysis of short-read metagenomic data using sourmash and YACHT.

## Acknowledgements

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE). We gratefully acknowledge the nf-core community for their tooling, framework, and support, which have been essential to the development of this pipeline.

> The nf-core framework for community-curated bioinformatics pipelines.  
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.  
> *Nat Biotechnol.* 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).

> Empowering bioinformatics communities with Nextflow and nf-core.  
> Langer, B.E., Amaral, A., Baudement, M.O., Bonath, F. et al.  
> *Genome Biol.* 26, 228 (2025). doi: [10.1186/s13059-025-03673-9](https://doi.org/10.1186/s13059-025-03673-9)