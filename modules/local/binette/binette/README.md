# binette/binette

Binette is a fast and accurate binning refinement tool designed to construct high-quality MAGs from the output of multiple binning tools.

## Description

Binette takes the output from multiple binning tools (e.g., MetaBAT2, MaxBin2, CONCOCT) and applies fundamental set operations (intersection, difference, union) to generate new hybrid bins. It then evaluates bin quality using CheckM2 to select the best possible bins.

Key features:
- **Enhanced Speed**: Runs CheckM2's initial steps once for all contigs, reusing intermediate results
- **No Limit on Input Bin Sets**: Supports any number of input bin sets, unlike metaWRAP
- **Dual Input Modes**: Supports both contig2bin tables and bin directories, similar to DAS Tool

## Usage

The module supports two input modes controlled by the `meta.input_mode` parameter:

### Mode 1: Contig2bin Tables

Use `meta.input_mode = 'contig2bin_tables'` when providing TSV files with contig-to-bin mappings:

```nextflow
include { BINETTE_BINETTE } from './modules/nf-core/binette/binette/main'

workflow {
    def meta = [ id:'sample1', input_mode:'contig2bin_tables' ]
    def contigs = file('assembly.fasta')
    def bins = [
        file('metabat2_bins.tsv'),
        file('maxbin2_bins.tsv')
    ]
    def proteins = [] // optional
    def checkm2_db = [] // optional
    
    BINETTE_BINETTE(
        tuple(meta, contigs, bins, proteins),
        checkm2_db
    )
}
```

### Mode 2: Bin Directories

Use `meta.input_mode = 'bin_dirs'` when providing directories containing bin FASTA files:

```nextflow
include { BINETTE_BINETTE } from './modules/nf-core/binette/binette/main'

workflow {
    def meta = [ id:'sample1', input_mode:'bin_dirs' ]
    def contigs = file('assembly.fasta')
    def bins = [
        file('metabat2_bins/', type: 'dir'),
        file('maxbin2_bins/', type: 'dir')
    ]
    def proteins = [] // optional
    def checkm2_db = [] // optional
    
    BINETTE_BINETTE(
        tuple(meta, contigs, bins, proteins),
        checkm2_db
    )
}
```

## Input

- `meta` (map): Groovy map containing sample information. Must include `input_mode` field set to either `'contig2bin_tables'` or `'bin_dirs'`
- `contigs` (file): Assembly contigs in FASTA format
- `bins` (files): Input bins in format determined by `meta.input_mode`:
  - If `contig2bin_tables`: List of TSV files with two columns (contig, bin)
  - If `bin_dirs`: List of directories containing bin FASTA files
- `proteins` (file, optional): Precomputed protein sequences in Prodigal format
- `checkm2_db` (file, optional): Path to CheckM2 diamond database

## Output

- `quality_reports` (file): TSV file with quality metrics for final bins
- `bins` (directory, optional): Directory containing final bins as FASTA files
- `contig2bin` (file): TSV file mapping contigs to final bin assignments
- `input_reports` (directory): Quality reports for input bin sets
- `log` (file): Binette execution log
- `versions` (file): Software versions

## Notes

- Binette requires a CheckM2 database. If not provided via `checkm2_db`, it will use the database configured with `checkm2 database`
- The `contamination_weight` parameter can be adjusted via `task.ext.args` to control the bin scoring function: `completeness - weight * contamination`
- Use `--no-write-fasta-bins` in `task.ext.args` to skip writing FASTA files for final bins

## References

- Mainguy, J., & Hoede, C. (2024). Binette: a fast and accurate bin refinement tool to construct high-quality Metagenome Assembled Genomes. *Journal of Open Source Software*, 9(102), 6782. https://doi.org/10.21105/joss.06782
- Chklovski, A., et al. (2023). CheckM2: a rapid, scalable, and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods*, 20(8), 1203-1212. https://doi.org/10.1038/s41592-023-01940-w