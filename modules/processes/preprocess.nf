#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_size_in_gb; swap_configurations; remove_alg_suffixes; add_suffix } from '../functions.nf'

process extract_b0 {
    memory { 4f * get_size_in_gb(dwi) }
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_b0") }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}__b0.nii.gz"), emit: b0
        tuple val(sid), path("${dwi.simpleName}__b0*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey b0 extract --in $dwi --bvals $bval --out ${dwi.simpleName}__b0 --config $config
        """
}

process squash_b0 {
    memory { 4f * get_size_in_gb(dwi) }
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}__b0_squashed.nii.gz"), path("${dwi.simpleName}__b0_squashed.bval"), path("${dwi.simpleName}__b0_squashed.bvec"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}__b0_squashed_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey b0 squash --in $dwi --bvals $bval --bvecs $bvec --out ${dwi.simpleName}__b0_squashed --config $config
        """
}
