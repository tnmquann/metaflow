#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { POSTPROCESS_READBASED } from '../main'

workflow {
    def postprocess_input_ch = channel.of([
        [id: 'default_test'],
        'default_test',
        'all',
        false,
        false,
        '--min_coverage 0.05 --use_average',
        file("${projectDir}/data/merged_sourmash_yacht.csv"),
        []
    ])
    def postprocess_script = file("${projectDir}/../../../../../bin/py_scripts/read_based/post_processing.py", checkIfExists: true)

    POSTPROCESS_READBASED(postprocess_input_ch, postprocess_script)

    POSTPROCESS_READBASED.out.default_output.view { _meta, default_dir ->
        "DEFAULT_DIR=${default_dir}"
    }
}
