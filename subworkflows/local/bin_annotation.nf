#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules for bin annotation
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA } from '../../modules/nf-core/bakta/bakta/main'
include { PROKKA } from '../../modules/nf-core/prokka/main'