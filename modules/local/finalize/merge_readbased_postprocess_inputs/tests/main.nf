#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MERGE_READBASED_POSTPROCESS_INPUTS } from '../main'

workflow {
    def aggregate_inputs = channel.of([
        [id: 'aggregate_test'],
        [
            file("${projectDir}/data/sample_a.csv"),
            file("${projectDir}/data/sample_b.csv")
        ]
    ])

    MERGE_READBASED_POSTPROCESS_INPUTS(aggregate_inputs)

MERGE_READBASED_POSTPROCESS_INPUTS.out.merged_csv.view { _meta, merged_csv ->
        "MERGED_CSV=${merged_csv}"
    }
}
