

process dmriqc_flow_input_tree {
    input:
        tuple val(sid), path(versa_input), path(versa_output)
    script:
        """
        # Produce Segmentation input
        mkdir -p Segment_Tissues
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/crop_wm_mask/${sid}_wm_mask__cropped.nii.gz \
            Segment_tissues/${sid}_mask_wm.nii.gz
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/crop_gm_mask/${sid}_gm_mask__cropped.nii.gz \
            Segment_tissues/${sid}_mask_gm.nii.gz
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/crop_csf_mask/${sid}_csf_mask__cropped.nii.gz \
            Segment_tissues/${sid}_mask_csf.nii.gz

        # Produce Bet DWI input
        mkdir -p Bet_DWI
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/t1_registration_wkf/transform_mask_to_b0/${sid}_t1_mask__resampled__transformed__transformed__transformed.nii.gz \
            Bet_DWI/${sid}_b0_bet_mask.nii.gz
        mkdir -p N4_DWI
        ln. -s $versa_output/all/$sid/preprocess/preprocess_wkf/resample_dwi/${sid}_dwi__to_eddy__eddy_corrected__checked__n4denoised__resampled.nii.gz \
            N4_DWI/${sid}_dwi_n4.nii.gz

        # Produce Bet T1 input
        mkdir -p Bet_T1
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/t1_registration_wkf/transform_mask_to_b0/${sid}_t1_mask__resampled__transformed__transformed__transformed.nii.gz \
            Bet_T1/${sid}_t1_bet_mask.nii.gz
        mkdir -p Resample_T1
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/t1_registration_wkf/transform_t1_to_b0/${sid}_t1__nlmeans_denoised__n4denoised__resampled__transformed.nii.gz \
            Bet_T1/${sid}_t1_resampled.nii.gz

        # Produce DWI denoise input
        mkdir -p Denoise_DWI
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/dwi_denoise_wkf/dwi_denoise/${sid}_dwi__even_dims__dwidenoised.nii.gz \
            Denoise_DWI/${sid}_dwi_denoised.nii.gz

        # Produce T1 denoise inputs
        mkdir -p Denoise_T1
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/t1_preprocess_wkf/nlmeans_denoise/${sid}_t1__nlmeans_denoised.nii.gz \
            Denoise_T1/${sid}_t1_denoised.nii.gz

        # Produce Gibbs correction input
        mkdir -p Gibbs_correction
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/dwi_gibbs_removal/${sid}_dwi__even_dims__dwidenoised__gibbs_corrected.nii.gz \
            Gibbs_correction/${sid}_dwi_gibbs_corrected.nii.gz
        mkdir -p Extract_DTI_Shell
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/check_odd_dimensions/${sid}_dwi__even_dims.bval \
            Extract_DTI_Shell/${sid}.bval
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/check_odd_dimensions/${sid}_dwi__even_dims.bvec \
            Extract_DTI_Shell/${sid}.bvec

        # Produce Eddy/Topup input
        mkdir -p Eddy_Topup
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/eddy_wkf/eddy/${sid}_dwi__to_eddy__eddy_corrected.bvec \
            Eddy_topup/${sid}.bvec
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/eddy_wkf/eddy/${sid}_dwi__to_eddy__eddy_corrected.nii.gz \
            Eddy_topup/${sid}_corrected.nii.gz
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/eddy_wkf/eddy/${sid}_dwi__to_eddy__eddy_corrected.bval \
            Eddy_topup/${sid}.bval_eddy
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/t1_mask_convert_datatype/${sid}_t1_mask__transformed__uint8.nii.gz \
            Eddy_Topup/${sid}_mask.nii.gz

        # Produce DWI resample input
        mkdir -p Resample_B0
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/extract_b0_preprocessed/${sid}_dwi__to_eddy__eddy_corrected__checked__n4denoised__resampled__cropped_b0.nii.gz \
            Resample_B0/${sid}_b0_resampled.nii.gz

        # Produce T1 resample input
        mkdir -p Resample_T1
        ln -s $versa_output/all/$sid/preprocess/preprocess_wkf/resample_t1/${sid}_t1__nlmeans_denoised__n4denoised__resampled.nii.gz \
            Resample_T1/${sid}_t1_resampled.nii.gz

        # Produce DTI metrics input
        mkdir -p DTI_Metrics
        ln -s $versa_output/$sid/dti/${sid}_dti_fa.nii.gz DTI_Metrics/${sid}_fa.nii.gz
        ln -s $versa_output/$sid/dti/${sid}_dti_md.nii.gz DTI_Metrics/${sid}_md.nii.gz
        ln -s $versa_output/$sid/dti/${sid}_dti_rd.nii.gz DTI_Metrics/${sid}_rd.nii.gz
        ln -s $versa_output/$sid/dti/${sid}_dti_ad.nii.gz DTI_Metrics/${sid}_ad.nii.gz
        ln -s $versa_output/$sid/dti/${sid}_dti_residuals.nii.gz DTI_Metrics/${sid}_residual.nii.gz
        ln -s $versa_output/$sid/dti/${sid}_dti_evecs_v1.nii.gz DTI_Metrics/${sid}_evecs_v1.nii.gz

        # Produce FRF input
        mkdir -p Compute_FRF
        ln -s $versa_output/all/$sid/reconstruct/reconstruct_wkf/csd_wkf/scilpy_response/${sid}_response.txt \
            Compute_FRF/${sid}_frf.txt

        # Produce FODF input
        mkdir -p FODF_Metrics
        scil_image_math.py addition \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_gm_afd.nii.gz \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_wm_afd.nii.gz
            FODF_Metrics/${sid}_afd_max.nii.gz -f
        scil_image_math.py addition \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_gm_afds.nii.gz \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_wm_afds.nii.gz
            FODF_Metrics/${sid}_afd_sum.nii.gz -f
        scil_image_math.py addition \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_gm_afdt.nii.gz \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_wm_afdt.nii.gz
            FODF_Metrics/${sid}_afd_total.nii.gz -f
        scil_image_math.py addition \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_gm_nufo.nii.gz \
            $versa_output/$sid/fodf/${sid}_fodf_metrics_wm_nufo.nii.gz
            FODF_Metrics/${sid}_nufo.nii.gz -f

        # Produce tracking input
        mkdir -p All_Tracking
        for tractogram in $versa_output/$sid/tracking/*.trk
        do
            ln -s \$tractogram All_Tracking/.
        done

        # Produce T1 registration input
        mkdir -p Register_T1
        ln -s $versa_output/$sid/${sid}_t1.nii.gz Register_T1/${sid}_t1_warped.nii.gz

        mkdir -p DTI_Metrics
        ln -s $versa_output/$sid/dti/${sid}_dti_rgb.nii.gz DTI_Metrics/${sid}_rgb.nii.gz

        # Produce segmentation input
        mkdir -p Segment_Tissues/
        ln -s $versa_output/$sid/segmentation/${sid}_wm_pvf.nii.gz Segment_Tissues/${sid}_map_wm.nii.gz
        ln -s $versa_output/$sid/segmentation/${sid}_gm_pvf.nii.gz Segment_Tissues/${sid}_map_gm.nii.gz
        ln -s $versa_output/$sid/segmentation/${sid}_csf_pvf.nii.gz Segment_Tissues/${sid}_map_csf.nii.gz

        mkdir -p Segment_Freesurfer/
        ln -s $versa_output/$sid/segmentation/${sid}_wm_mask.nii.gz Segment_Freesurfer/${sid}_mask_wm.nii.gz
        ln -s $versa_output/$sid/segmentation/${sid}_gm_mask.nii.gz Segment_Freesurfer/${sid}_mask_gm.nii.gz
        ln -s $versa_output/$sid/segmentation/${sid}_csf_mask.nii.gz Segment_Freesurfer/${sid}_mask_csf.nii.gz

        # Produce PFT maps input
        mkdir -p PFT_Seeding_Mask
        ln -s $versa_output/$sid/tracking/${sid}_wm_gm_interface.nii.gz PFT_Seeding_Mask/${sid}_seeding_mask.nii.gz

        mkdir -p PFT_Tracking_Maps
        ln -s $versa_output/$sid/tracking/${sid}_map_include.nii.gz PFT_Tracking_Maps/${sid}_map_include.nii.gz
        ln -s $versa_output/$sid/tracking/${sid}_map_exclude.nii.gz PFT_Tracking_Maps/${sid}_map_exclude.nii.gz

        # Produce input QC tree
        mkdir -p Inputs_QC
        ln -s $versa_input/$sid/${sid}_dwi.bval Inputs_QC/$sid/dwi.bval
        [ -f $versa_input/$sid/${sid}_rev.bval ] && ln -s $versa_input/$sid/${sid}_rev.bval Inputs_QC/$sid/rev.bval
        ln -s $versa_input/$sid/${sid}_dwi.bvec Inputs_QC/$sid/dwi.bvec
        [ -f $versa_input/$sid/${sid}_rev.bvec ] && ln -s $versa_input/$sid/${sid}_rev.bval Inputs_QC/$sid/rev.bvec
        ln -s $versa_input/$sid/${sid}_t1.nii.gz Inputs_QC/$sid/t1.nii.gz
        ln -s $versa_input/$sid/${sid}_dwi.nii.gz Inputs_QC/$sid/dwi.nii.gz
        [ -f $versa_input/$sid/${sid}_rev.bval ] && ln -s $versa_input/$sid/${sid}_rev.nii.gz Inputs_QC/$sid/rev_dwi.nii.gz
        """
}