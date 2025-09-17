#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { METAFLOW } from './workflow/metaflow'

workflow {
    METAFLOW()
}