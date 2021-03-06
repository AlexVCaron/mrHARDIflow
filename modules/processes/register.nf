#!/usr/bin/env nextflow

import java.io.File

nextflow.enable.dsl=2

params.random_seed = 1234


include { remove_alg_suffixes } from '../functions.nf'

process ants_register {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : f.contains("registration_warped.nii.gz") ? publish_suffix ? "${sid}_${publish_suffix}.nii.gz" : remove_alg_suffixes(f) : null : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(moving), path(target), val(reference), file(mask), file(metadata)
        val(caller_name)
        val(additional_publish_path)
        val(publish)
        val(publish_suffix)
        path(config)
    output:
        tuple val(sid), path("${moving[0].simpleName}__registration_ref.nii.gz"), emit: reference
        tuple val(sid), path("${moving[0].simpleName}__[A-Z]_registration_*.*"), emit: transformation
        tuple val(sid), path("${moving[0].simpleName}__[A-Z]_inv_registration_*.*"), emit: inverse_transformation
        tuple val(sid), path("${moving[0].simpleName}__registration_warped.nii.gz"), optional: true, emit: image
        tuple val(sid), path("${moving[0].simpleName}__registration_warped_metadata.*"), optional: true, emit: metadata
    script:
        def mask_arg = ""
        if ( !mask.iterator().inject(false) { c, i -> c || i.empty() } ) {
            mask_arg = "--mask ${mask.iterator().collect{ it.name }.join(',')}"
        }

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed

        magic-monkey ants_registration \
            --moving ${moving.join(",")} \
            --target ${target.join(",")} \
            --out ${moving[0].simpleName}__registration $mask_arg \
            ${params.verbose_outputs ? "--verbose" : ""} \
            --config $config

        cp ${file(reference).name} ${moving[0].simpleName}__registration_ref.nii.gz
        cnt1=0
        cnt2=1
        while true
        do
            found=false
            warp_found=false
            if [ -f ${moving[0].simpleName}__registration\${cnt1}GenericRigid.mat ]
            then
                printf -v letter "\\x\$(printf %x \$((\$cnt2 + 64)))"
                (( ++cnt2 ))
                cp ${moving[0].simpleName}__registration\${cnt1}GenericRigid.mat ${moving[0].simpleName}__\${letter}_registration_rigid.mat
                cp ${moving[0].simpleName}__registration\${cnt1}GenericRigid.mat ${moving[0].simpleName}__\${letter}_inv_registration_rigid.mat
                found=true
            fi
            if [ -f ${moving[0].simpleName}__registration\${cnt1}GenericAffine.mat ]
            then
                printf -v letter "\\x\$(printf %x \$((\$cnt2 + 64)))"
                (( ++cnt2 ))
                cp ${moving[0].simpleName}__registration\${cnt1}GenericAffine.mat ${moving[0].simpleName}__\${letter}_registration_affine.mat
                cp ${moving[0].simpleName}__registration\${cnt1}GenericAffine.mat ${moving[0].simpleName}__\${letter}_inv_registration_affine.mat
                found=true
            fi
            if [ -f ${moving[0].simpleName}__registration\${cnt1}Warp.nii.gz ]
            then
                printf -v letter "\\x\$(printf %x \$((\$cnt2 + 64)))"
                (( ++cnt2 ))
                mv ${moving[0].simpleName}__registration\${cnt1}Warp.nii.gz ${moving[0].simpleName}__\${letter}_registration_syn.nii.gz
                mv ${moving[0].simpleName}__registration\${cnt1}InverseWarp.nii.gz ${moving[0].simpleName}__\${letter}_inv_registration_syn_inverse.nii.gz
                found=true
            fi
            
            if \$found
            then
                (( ++cnt1 ))
            else
                break
            fi
        done
            
        """
}

process ants_correct_motion {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(moving), path(target), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${moving[0].simpleName}__motion_correct_warped.nii.gz"), emit: image
        tuple val(sid), path("${moving[0].simpleName}__motion_correct_warped_metadata.*"), optional: true, emit: metadata
    script:
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed
        magic-monkey ants_motion --moving ${moving.join(",")} --target ${target.join(",")} --out ${sid}__motion_correct --config $config
        """
}

process ants_transform {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : publish_suffix ? "${sid}_${publish_suffix}.nii.gz" : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(img), path(ref), path(trans), val(invert), file(bvec), file(metadata)
        val(caller_name)
        val(additional_publish_path)
        val(publish)
        val(publish_suffix)
        path(config)
    output:
        tuple val(sid), path("${img.simpleName}__transformed.nii.gz"), emit: image
        tuple val(sid), path("${img.simpleName}__transformed.bvec"), optional: true, emit: bvec
        tuple val(sid), path("${img.simpleName}__transformed_metadata.*"), optional: true, emit: metadata
    script:
        def args = "--in $img --ref $ref"
        def trans_str = (trans instanceof Path) ? trans : trans.join(',')
        if ( trans && (trans instanceof Path) ? !trans.empty() : !trans.isEmpty() ) {
            args += " --trans $trans_str"
        }
        if ( !bvec.empty() ) {
            args += " --bvecs $bvec"
        }
        if ("$invert") {
            args += " --inv ${invert.join(",")}"
        }
        """
        export ANTS_RANDOM_SEED=$params.random_seed
        magic-monkey ants_transform $args --out ${img.simpleName}__transformed --config $config
        """
}
