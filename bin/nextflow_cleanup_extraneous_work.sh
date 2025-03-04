#! /bin/bash

## See also https://github.com/nextflow-io/nextflow/discussions/4308

## cd to a parent directory for a Nextflow pipeline executation, i.e. contains .nextflow and work directories

## Find work directories essential to the last pipeline run, as absolute paths
nextflow log last > /tmp/preserve_dirs.txt

## Find all work directories, as absolute paths
find "$(readlink -f ./work)" -maxdepth 2 -type d -path '**/work/*/*' > /tmp/all_dirs.txt

## Concatenate, sort, and count, filtering to those that show up only once (i.e., just once from all_dirs.txt)
cat /tmp/all_dirs.txt /tmp/preserve_dirs.txt | sort | uniq -c | grep "      1 /" | grep -Po '/.+' > /tmp/to_delete_dirs.txt

## Delete the extraneous work directories
cat /tmp/to_delete_dirs.txt | xargs -r -P 4 -n 1 rm -rf

## Clean up
rm -f /tmp/preserve_dirs.txt /tmp/all_dirs.txt /tmp/to_delete_dirs.txt
