#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false

include { load_dataset } from "./workflows/io.nf"
include { register_atlases_wkf } from "./workflows/register.nf"

params.data_input = false

def get_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}

workflow {

    t1_channel = Channel.fromFilePairs(
        "${params.data_input}/**/*_t1.nii.gz", maxDepth: 2, size: 1
    ){ get_id(it.parent, root) }
    brainmask_channel = Channel.fromFilePairs(
        "${params.data_input}/**/*_brainmask.nii.gz", maxDepth: 2, size: 1
    ){ get_id(it.parent, root) }
    gm_mask_channel = Channel.fromFilePairs(
        "${params.data_input}/**/segmentation/*_gm_mask.nii.gz", maxDepth: 3, size: 1
    ){ get_id(it.parent, root) }
    transforms_channel = Channel.fromFilePairs(
        "${params.data_input}/**/transforms/b0_to_template/scripts_transforms/*image*b0_to_template*.nii.gz"
    ){ get_id(it.parent, root) }

    register_atlases_wkf(t1_channel, brainmask_channel, gm_mask_channel, transforms_channel)
}

workflow.onComplete {
    log.info "Pipeline completed at : $workflow.complete"
    log.info "Execution status : ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration : $workflow.duration"
}
