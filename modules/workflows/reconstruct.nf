#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { diamond; mrtrix_dti; csd; response; scilpy_response; scilpy_msmt_response; scilpy_csd; scilpy_msmt_csd } from '../processes/reconstruct.nf'
include { scil_dti_and_metrics } from '../processes/measure.nf'
include { tournier2descoteaux_odf; extract_shells; check_for_duplicates } from '../processes/utils.nf'
include { get_config_path } from '../functions.nf'

params.reconstruct_use_mrtrix = false
params.convert_tournier2descoteaux = true
params.msmt_odf = false

params.reconstruct_diamond_config = file("${get_config_path()}/reconstruct_diamond_config.py")
params.reconstruct_mrtrix_dti_config = file("${get_config_path()}/reconstruct_mrtrix_dti_config.py")
params.reconstruct_mrtrix_csd_config = file("${get_config_path()}/reconstruct_mrtrix_csd_config.py")
params.reconstruct_mrtrix_frf_config = file("${get_config_path()}/reconstruct_mrtrix_frf_config.py")
params.extract_shell_greater_than_one_config = file("${get_config_path()}/extract_shell_greater_than_one_config.py")


workflow csd_wkf {
    take:
        dwi_channel
        mask_channel
        tissue_masks_channel
        safe_wm_mask_channel
    main:
        response_channel = Channel.empty()
        odfs_channel = Channel.empty()

        check_for_duplicates(dwi_channel.map{ it + [""] }, "reconstruct")
        dwi_channel = check_for_duplicates.out.dwi

        if ( params.reconstruct_use_mrtrix && !params.msmt_odf ) {
            response(dwi_channel.join(mask_channel), "reconstruct", params.reconstruct_mrtrix_frf_config)
            csd(response.out.responses.join(dwi_channel.join(mask_channel)), "reconstruct", params.reconstruct_mrtrix_csd_config)
            response_channel = response.out.responses
            odfs_channel = csd.out.odfs
            if ( params.convert_tournier2descoteaux ) {
                tournier2descoteaux_odf(csd.out.odfs, "reconstruct")
                csd_channel = tournier2descoteaux_odf.out.odfs
            }
        }
        else {
            if ( params.msmt_odf ) {
                dwi_channel = extract_shells(dwi_channel, "reconstruct", params.extract_shell_greater_than_one_config)
                scilpy_msmt_response(dwi_channel.join(mask_channel).join(tissue_masks_channel), "reconstruct")
                scilpy_msmt_csd(dwi_channel.join(scilpy_msmt_response.out.response).join(mask_channel), "reconstruct")
                response_channel = scilpy_msmt_response.out.response
                odfs_channel = scilpy_msmt_csd.out.odfs
            }
            else {
                scilpy_response(dwi_channel.join(mask_channel).join(safe_wm_mask_channel), "reconstruct")
                scilpy_csd(dwi_channel.join(scilpy_response.out.response).join(mask_channel), "reconstruct")
                response_channel = scilpy_response.out.response
                odfs_channel = scilpy_csd.out.odfs
            }
        }
    emit:
        odfs = odfs_channel
        responses = response_channel
}

workflow dti_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        dti_output = Channel.empty()
        if ( params.reconstruct_use_mrtrix ) {
            mrtrix_dti(dwi_channel.join(mask_channel), "reconstruct", params.reconstruct_mrtrix_dti_config)
            dti_output = mrtrix_dti.out.dti
        }
        else {
            scil_dti_and_metrics(dwi_channel.join(mask_channel), "reconstruct", "measure")
            dti_output = scil_dti_and_metrics.out.dti
        }
    emit:
        dti = dti_output
}

workflow diamond_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        dwi_channel = dwi_channel.groupTuple()
        dwi_image = dwi_channel.map{ [it[0], it[1]] }
        other_files = dwi_channel.map{ [it[0], it.subList(2, it.size()).inject([]){ c, t -> c + t }] }
        diamond(dwi_image.join(mask_channel).join(other_files), "reconstruct", params.reconstruct_diamond_config)
    emit:
        data = diamond.out.diamond
        xml_summary = diamond.out.xml_summary
}
