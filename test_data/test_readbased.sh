#! /bin/bash

# Please check ./conf/test.config for customization options.
projectDir="/path/to/your/main_nf" && \
export projectDir && \
envsubst < ${projectDir}/test_data/sample_e2e.csv > /tmp/sample_temp.csv && \
envsubst < ${projectDir}/test_data/mockupdb_r226_k31_thresh_0.995/mockupdb_r226_k31_thresh_0.995_config.json > /tmp/config_temp.json && \
nextflow run ${projectDir}/main.nf -profile test,conda --output path/to/output