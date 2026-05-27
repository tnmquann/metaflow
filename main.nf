#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { METAFLOW } from './workflow/metaflow'

workflow {
    if (!(workflow.preview && !params.input && !params.assembly_input)) {
        METAFLOW()
    }
}
