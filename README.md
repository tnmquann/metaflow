# profiler_sourmash
A read-based shotgun metagenomics processing pipeline built on the Nextflow DSL2 framework.

## Parameters
| Parameter                 | Description                                    | Default value                                                               |
| ------------------------- | ---------------------------------------------- | --------------------------------------------------------------------------- |
| `base_dir`                | Base directory for workflow data               | ` `                                                                         |
| `workDir`                 | Working directory for Nextflow                 | ` `                                                                         |
| `trimmed_fastq`           | Input directory containing trimmed fastq files | ` `                                                                         |
| `outdir`                  | Output directory for results                   | `${params.base_dir}/results`                                                |
| `ksize`                   | K-mer size for sourmash                        | `51`                                                                        |
| `sourmash_database`       | Path to sourmash database                      | `Please use RocksDB format `                                                |
| `yacht_database`          | Path to yacht database                         | `path/to/json/yacht/file`                                                   |
| `enable_copymergedseqs`   | Enable copying of merged sequences             | `false`                                                                     |
| `enable_copysketch`       | Enable copying of sketch files                 | `true`                                                                      |
| `enable_copyintermediate` | Enable copying of intermediate files           | `true`                                                                      |
| `cleanup`                 | Enable cleanup of temporary files              | `false`                                                                     |
| `cleanup_patterns`        | File patterns to clean up                      | `['work', 'intermediate_files', '.nextflow*', '.command.*', 'process_tmp']` |
| `publish_dir_mode`        | Mode for publishing output files               | `copy`                                                                      |
| `email`                   | Email for notifications                        | `null`                                                                      |
| `cleanup_dry_run`         | Preview cleanup without deleting               | `false`                                                                     |
| `keep_days`               | Keep files newer than specified days           | `7`                                                                         |

## How to run
```
nextflow run main.nf -profile conda
```