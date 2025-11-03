#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules for bin taxonomic classification
include { GTDBTK_CLASSIFYWF } from '../../modules/nf-core/gtdbtk/classifywf/main'

// Import local modules for bin taxonomic classification
include { GTDBTK_DB_PREPARATION } from '../../modules/local/gtdbtk/dbpreparation/main'
include { GTDBTK_SUMMARY } from '../../modules/local/gtdbtk/summary/main'