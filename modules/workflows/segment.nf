#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    registration_wkf as nmt_registration_wkf;
    registration_wkf as wm_seg_registration_wkf
} from "./preprocess.nf"
include { propagate_labels_wkf as propagate_atlas } from './labels.nf'
include { atropos } from '../processes/segment.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'
include { 
    prepend_sid as prepend_sid_template;
    prepend_sid as prepend_sid_template_mask;
    prepend_sid as prepend_sid_segmentation;
    prepend_sid as prepend_sid_template_fa;
    prepend_sid as prepend_sid_wm_atlas;
    prepend_sid as prepend_sid_d99;
    prepend_sid as prepend_sid_charm;
    prepend_sid as prepend_sid_sarm;
    prepend_sid as prepend_sid_inia19
    pvf_to_mask
} from '../processes/utils.nf'
include {
    resampling_reference;
    resampling_reference as resampling_reference_fa;
    scilpy_resample_to_reference as resample_template;
    scilpy_resample_to_reference as resample_segmentation;
    scilpy_resample_to_reference as resample_segmentation_pre_transform;
    scilpy_resample_to_reference as resample_template_fa;
    scilpy_resample_to_reference as resample_wm_atlas;
    scilpy_resample_to_reference as resample_atlas;
} from '../processes/upsample.nf'
include {
    ants_transform as transform_atlas;
    ants_transform as transform_nmt
} from '../processes/register.nf'
include { get_config_path; get_data_path } from '../functions.nf'   

params.force_resampling_resolution = false
params.resampling_min_resolution = false
params.resampling_subdivision = 2

params.segmentation_classes = ["csf", "wm", "gm", "dgm", "pdgm", "blood"]
params.atropos_n4_tissues = ["wm", "dgm", "pdgm"]
params.csf_extent_lambda = 0.1
params.csf_extent_probability = 0.5
params.default_extent_lambda = 0.1
params.default_extent_probability = 0.8
params.default_blur_factor = 0.4
params.csf_blur_factor = 0.25

params.tissue_segmentation_root = "${get_data_path()}/maccaca_mulatta/tissue_segmentation"
params.wm_segmentation_root = "${get_data_path()}/maccaca_mulatta/wm_segmentation"

params.segmentation_registration_config = file("${get_config_path()}/segmentation_registration_config.py")
params.ants_transform_segmentation_config = file("${get_config_path()}/ants_transform_segmentation_config.py")


workflow segment_nmt_wkf {
    take:
        t1_channel
        mask_channel
        template_resampling_reference
        template_to_t1_transform
    main:
        segmentation_channel = prepend_sid_segmentation(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation.nii.gz")] })

        if ( template_to_t1_transform && template_resampling_reference ) {
            resample_segmentation_pre_transform(
                segmentation_channel
                    .join(template_resampling_reference)
                    .map{ it + ["", ""] },
                "segmentation",
                "nn",
                false,
                false,
                "", ""
            )

            transform_nmt(
                resample_segmentation_pre_transform.out.image
                    .join(t1_channel)
                    .join(template_to_t1_transform)
                    .map{ it + ["", ""] },
                "segmentation",
                "", "false", "",
                params.ants_transform_segmentation_config
            )

            nmt_in_subject_space = transform_nmt.out.image
        }
        else {
            template_channel = prepend_sid_template(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_t1.nii.gz")] })
            template_mask_channel = prepend_sid_template_mask(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_mask_whole_no_bv.nii.gz")] })

            resampling_reference(
                t1_channel.join(template_channel).map{ [it[0], it[1..-1]] },
                "segmentation",
                params.resampling_subdivision,
                params.resampling_min_resolution,
                ""
            )
            resample_template(
                template_channel
                    .join(resampling_reference.out.reference)
                    .join(template_mask_channel)
                    .map{ it + [""] },
                "segmentation",
                "lin",
                false,
                false,
                "", ""
            )
            resample_segmentation(
                segmentation_channel
                    .join(resampling_reference.out.reference)
                    .map{ it + ["", ""] },
                "segmentation",
                "nn",
                false,
                false,
                "", ""
            )

            nmt_registration_wkf(
                t1_channel.map{ [it[0], [it[1]]] },
                resample_template.out.image.map{ [it[0], [it[1]]] },
                resample_segmentation.out.image.map{ [it[0], [it[1]]] },
                mask_channel.join(resample_template.out.mask).map{ [it[0], [it[1]]] },
                null,
                null,
                "segmentation",
                false,
                "", "",
                params.segmentation_registration_config,
                params.ants_transform_segmentation_config
            )

            template_resampling_reference = resampling_reference.out.reference
            template_to_t1_transform = nmt_registration_wkf.out.transform
            nmt_in_subject_space = nmt_registration_wkf.out.image
        }

        def class_mapping = params.segmentation_classes.inject([:]){map, class ->
            map[class] = [
                label: params.segmentation_classes.indexOf(class) + 1,
                blur: is_csf ? params.csf_blur_factor : params.default_blur_factor,
                lambda: is_csf ? params.csf_extent_lambda : params.default_extent_lambda,
                pdist: is_csf ? params.csf_extent_probability : params.default_extent_probability,
                n4: params.atropos_n4_tissues.contains(class)
            ]
        }

        atropos(
            t1_channel.join(mask_channel).join(nmt_in_subject_space),
            class_mapping,
            "segmentation"
        )

        pvf_to_mask(
            atropos.out.vol_fractions
                .map{ [it[0], it[1].sort{ a, b -> params.segmentation_classes.findIndexOf{ i -> a.simpleName.contains("_$i") } <=> params.segmentation_classes.findIndexOf{ i -> b.simpleName.contains("_$i") } }] }
                .join(mask_channel),
            "segmentation",
            "segmentation"
        )

    emit:
        segmentation = atropos.out.segmentation
        volume_fractions_full = atropos.out.vol_fractions
        volume_fractions = pvf_to_mask.out.pvf_3t
        tissue_masks = pvf_to_mask.out.masks
            .map{ [ it[0], it[1].sort{ a, b -> ["wm", "gm", "csf"].findIndexOf{ i -> a.simpleName.contains(i) } <=> ["wm", "gm", "csf"].findIndexOf{ i -> b.simpleName.contains(i) } }] }
        safe_wm_mask = pvf_to_mask.out.safe_wm_mask
}

workflow segment_atlas_wkf {
    take:
        atlas_channel
        reference_channel
        transform_channel
        mask_channel
        resegment_atlas
        class_mapping
        publish_path
    main:
        transformed_atlas = Channel.empty()

        if ( resegment_atlas ) {
            propagate_atlas(
                reference_channel,
                mask_channel,
                atlas_channel,
                transform_channel,
                class_mapping,
                publish_path
            )

            transformed_atlas = propagate_atlas.out.labels
        }
        else {
            transform_atlas(
                atlas_channel
                    .join(reference_channel)
                    .join(transform_channel)
                    .map{ it + ["", ""] },
                "atlases", "atlases", "true", "",
                params.ants_transform_segmentation_config
            )

            transformed_atlas = transform_atlas.out.image
        }
    emit:
        atlas = transformed_atlas
}

workflow segment_wm_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        template_fa_channel = prepend_sid_template_fa(dwi_channel.map{ [it[0], file("${params.wm_segmentation_root}/wm_segmentation_fa.nii.gz")] })
        wm_atlas_channel = prepend_sid_wm_atlas(dwi_channel.map{ [it[0], file("${params.wm_segmentation_root}/wm_segmentation_atlas.nii.gz")] })

        resampling_reference_fa(
            dwi_channel.map{ it[0..1] }.join(template_fa_channel).map{ [it[0], it[1..-1]] },
            "segmentation",
            params.resampling_subdivision,
            params.resampling_min_resolution,
            ""
        )
        resample_template_fa(
            template_fa_channel
                .join(resampling_reference_fa.out.reference)
                .map{ it + ["", ""] },
            "preprocess",
            "lin",
            false,
            false,
            "", ""
        )
        resample_wm_atlas(
            wm_atlas_channel
                .join(resampling_reference_fa.out.reference)
                .map{ it + ["", ""] },
            "preprocess",
            "nn",
            false,
            false,
            "", ""
        )

        scil_compute_dti_fa(dwi_channel.join(mask_channel), "segmentation", "segmentation", false)

        wm_seg_registration_wkf(
            scil_compute_dti_fa.out.fa.map{ [it[0], [it[1]]] },
            resample_template_fa.out.image.map{ [it[0], [it[1]]] },
            resample_wm_atlas.out.image.map{ [it[0], [it[1]]] },
            null,
            null,
            null,
            "segmentation",
            true,
            "", "",
            params.segmentation_registration_config,
            params.ants_transform_segmentation_config
        )
    emit:
        segmentation = wm_seg_registration_wkf.out.image
}
