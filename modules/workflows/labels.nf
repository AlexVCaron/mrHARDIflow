#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { atropos } from '../processes/segment.nf'

include {ants_transform as transform_labels} from '../processes/register.nf'

workflow propagate_labels_wkf {
    take:
        t1_channel
        mask_channel
        labels_channel
        affine_transform_channel
        class_mapping
        publish_path
    main:
        transform_labels(
            labels_channel
                .join(t1_channel)
                .join(affine_transform_channel)
                .map{ it + ["", ""] },
            "", "", "false", "",
            params.ants_transform_segmentation_config
        )

        atropos(
            t1_channel
                .join(mask_channel)
                .join(transform_labels.out.image),
            class_mapping,
            publish_path
        )
    emit:
        labels = atropos.out.segmentation
}