#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.register_d99 = true
params.register_charm = true
params.register_sarm = true
params.register_inia19 = true
params.resegment_atlases = false

params.atlas_resegment_blur = 0.5
params.atlas_resegment_lambda = 0.1
params.atlas_resegment_extent_probability = 0.8

params.tissue_segmentation_root = "${get_data_path()}/maccaca_mulatta/tissue_segmentation"

workflow register_atlases_wkf {
    take:
        reference_channel
        transform_channel
        mask_channel
        gm_mask_channel
    main:
        transformed_d99 = Channel.empty()
        transformed_charm = Channel.empty()
        transformed_sarm = Channel.empty()
        transformed_inia19 = Channel.empty()

        if ( params.register_d99 ) {
            d99_channel = prepend_sid_d99(
                reference_channel
                    .map{ [it[0], file("${params.tissue_segmentation_root}/D99_atlas.nii.gz")] }
            )

            d99_class_mapping = file("${params.tissue_segmentation_root}/D99_labels.txt")
                .readLines()
                .inject([:]){ map, label ->
                map[label] = [
                    label: label,
                    blur: params.atlas_resegment_blur,
                    lambda: params.atlas_resegment_lambda,
                    pdist: params.atlas_resegment_extent_probability
                ]
            }

            segment_d99(
                d99_atlas,
                reference_channel,
                transform_channel,
                gm_mask_channel,
                params.resegment_atlases,
                d99_class_mapping,
                "atlases"
            )

            transformed_d99 = segment_d99.out.atlas
        }

        /*if ( params.register_charm ) {
            charm_channel = prepend_sid_charm(
                template_resampling_reference
                    .map{ [it[0], file("${params.tissue_segmentation_root}/CHARM_atlas.nii.gz")] }
            )
            
            segment_charm(
                charm_atlas,
                resampling_reference_channel,
                registration_reference_channel,
                registration_transform_channel,
                mask_channel,
                params.resegment_atlases,
                params.charm_class_mapping,
                "atlases"
            )

            transformed_charm = segment_charm.out.atlas
        }

        if ( params.register_sarm ) {
            sarm_channel = prepend_sid_sarm(
                template_resampling_reference
                    .map{ [it[0], file("${params.tissue_segmentation_root}/SARM_atlas.nii.gz")] }
            )
            
            segment_sarm(
                sarm_atlas,
                resampling_reference_channel,
                registration_reference_channel,
                registration_transform_channel,
                mask_channel,
                params.resegment_atlases,
                params.sarm_class_mapping,
                "atlases"
            )

            transformed_sarm = segment_sarm.out.atlas
        }*/

        if ( params.register_inia19 ) {
            inia19_channel = prepend_sid_inia19(
                reference_channel
                    .map{ [it[0], file("${params.tissue_segmentation_root}/INIA19_atlas.nii.gz")] }
            )

            inia19_class_mapping = file("${params.tissue_segmentation_root}/INIA19_labels.txt")
                .readLines()
                .inject([:]){ map, label ->
                map[label] = [
                    label: label,
                    blur: params.atlas_resegment_blur,
                    lambda: params.atlas_resegment_lambda,
                    pdist: params.atlas_resegment_extent_probability
                ]
            }

            segment_inia19(
                inia19_atlas,
                reference_channel,
                transform_channel,
                mask_channel,
                params.resegment_atlases,
                inia19_class_mapping,
                "atlases"
            )

            transformed_inia19 = segment_inia19.out.atlas
        }
    emit:
        d99_atlas = transformed_d99
        charm_atlas = transformed_charm
        sarm_atlas = transformed_sarm
        inia19_atlas = transformed_inia19
}