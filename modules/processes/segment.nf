#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.atropos_prior_weight = 0.2
params.atropos_prior_threshold = false
params.atropos_mrf_weight = 0.3
params.atropos_mrf_neighborhood = 1

params.atropos_niter = 5
params.atropos_convergence = 0.001
params.atropos_n4_niter = 15
params.atropos_outliers_window = "BoxPlot[0.25,0.75,1.5]"

params.atropos_adaptive_smoothing = false
params.atropos_bspline_levels = 6
params.atropos_bspline_order = 3
params.atropos_bspline_mesh_factor = 1.0
params.atropos_euclidean = false

params.atropos_prior = "PriorLabelImage"
params.atropos_posterior = "Aristotle"
params.atropos_likelihood = "Gaussian"
params.atropos_use_mixture_model = true
params.atropos_inital_temperature = 1.0
params.atropos_minimum_temperature = 0.1
params.atropos_annealing_rate = 1.0

params.atropos_n4_iterations = [200, 100, 50, 50]
params.atropos_n4_convergence_eps = 1E-10
params.atropos_n4_shrink_factor = 2
params.atropos_n4_bspline_fitting = 100

params.random_seed = 1234
params.disable_n4 = false

include { remove_alg_suffixes } from '../functions.nf'


/*process prepare_atropos {
    label "LIGHTSPEED"
    label "res_single_cpu"

    input:
        tuple val(sid), path(anat_images), path(segmentation)
    output:
        tuple val(sid), path("${sid}_run_atropos.sh"), emit: atropos_script
        tuple val(sid), path("${sid}_*_priors.nii.gz"), emit: segmentation_priors
    script:
        def args = ""
        if (!tissue_mappings.empty()) args += " --mappings $tissue_mappings"
        """
        mrhardi seg2pvf $args \
            --in $segmentation \
            --wm-label $params.segmentation_classes.indexOf("wm") \
            --gm-label $params.segmentation_classes.indexOf("gm") \
            
            --prefix ${segmentation.simpleName}_
        """
}*/


process atropos {
    label "ATROPOS"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/segmentation", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(t1_image), path(mask), path(segmentation)
        val(class_mapping) /*
            [ tissue_class: [ label: int, blur: float, lambda: float, pdist: float, n4: boolean ] ]
        */
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_segmentation.nii.gz"), emit: segmentation
        tuple val(sid), path("${sid}_*_pvf.nii.gz"), emit: vol_fractions
    script:
        class_mapping = class_mapping.sort{ it1, it2 -> it1.value.label <=> it2.value.label }
        def basename = segmentation.simpleName
        def classes = class_mapping.values()
        def tissues = class_mapping.keySet()
        def labels = class_mapping.collect{ it.value.label.toString() }
        def width = labels[-1].length()
        def n4_labels = class_mapping.findAll{ it.value.n4 }.collect{ it.value.label }

        def prior_init = ""
        if (params.disable_n4 && params.atropos_prior == "PriorLabelImage") {
            prior_init = [classes.size(), segmentation,$params.atropos_prior_weight].join(",")
        }
        else if (params.disable_n4 && params.atropos_prior == "PriorProbabilityImages") {
            prior_init = [
                classes.size(),
                "${basename}_%0${width}d.nii.gz",
                params.atropos_prior_weight,
                params.atropos_prior_threshold
            ].join(",")
        }

        def posterior = [
            params.atropos_use_mixture_model ? 1 : 0,
            params.atropos_inital_temperature,
            params.atropos_minimum_temperature,
            params.atropos_annealing_rate
        ].join(",")

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OMP_NUM_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed

        function join_array { local IFS="\$1"; shift; echo "\$*"; }

        tissues=( ${tissues.collect{ "\"$it\"" }.join(' ')} )
        values=( ${labels.collect{ "\"$it\"" }.join(' ')} )
        labels=( ${labels.collect{ "\"${ '0' * (width - it.length()) + it }\"" }.join(' ')} )
        blur=( ${class_mapping.collect{ "\"${ it.value.blur }\"" }.join(' ')} )

        mrhardi seg2mask \
            --in $segmentation \
            --values \$(join_array , \${values[@]}) \
            --labels \$(join_array , \${labels[@]}) \
            --out $basename

        if [ -f $mask ]
        then

        mask=$mask

        else

        python3 << - ENDSCRIPT
        import nibabel as nib
        img = nib.load("$segmentation")
        mask = (img.get_fdata() > 0).astype(np.uint8)
        nib.save(nib.Nifti1Image(mask, img.affine), "mask.nii.gz")
        ENDSCRIPT

        mask=mask.nii.gz

        fi

        for ix in \${!labels[@]}
        do

        scil_image_math.py blur ${basename}_\${labels[ix]}.nii.gz \
            \${blur[ix]} \
            ${basename}_\${labels[ix]}.nii.gz \
            --data_type float32 -f

        done

        spacing=\$(mrinfo -spacing $t1_image | awk '{print \$1}')

        if [ $params.disable_n4 = true ]
        then

        mesh_size=\$(echo "\${spacing[0]}*$params.atropos_bspline_mesh_factor" | bc)
        mesh_size=\$(printf "%.2fx%.2fx%.2f" \$mesh_size \$mesh_size \$mesh_size)

        Atropos -d 3 - r 0 -a [$t1_image,$params.atropos_adaptive_smoothing] -x \$mask \
            -b [$params.atropos_bspline_levels,$mesh_size,$params.atropos_bspline_order] \
            -i ${params.atropos_prior}[$prior_init] \
            -p ${params.atropos_posterior}[$posterior] \
            -c [${params.atropos_niter},${params.atropos_convergence}] \
            -k $params.atropos_likelihood \
            -m [$params.atropos_mrf_weight,${([params.atropos_mrf_neighborhood] * 3).join('x')})] \
            -o [atropos_segmentation.nii.gz,atropos_SegmentationPosteriors%0${width}d.nii.gz] \
            ${classes.collect{ "-l ${it.label}[${it.lambda},${it.pdist}]" }.join(' \\\n    ')}
            -w $params.atropos_outliers_window \
            -e ${params.atropos_euclidean ? 1 : 0}

        else

        antsAtroposN4.sh -u 0 -d 3 \
            -a $t1_image \
            -x \$mask \
            -m $params.atropos_n4_niter \
            -n $params.atropos_niter \
            -b ${params.atropos_posterior}[${params.atropos_use_mixture_model ? 1 : 0}] \
            ${classes.collect{ "-l ${it.label}[${it.lambda},${it.pdist}]" }.join(' \\\n    ')}
            -r [$params.atropos_mrf_weight,${([params.atropos_mrf_neighborhood] * 3).join('x')}] \
            ${n4_labels.collect{ "-y ${it}" }.join(' ')} \
            -e [${params.atropos_n4_iterations.join('x')},$params.atropos_n4_convergence_eps] \
            -f $params.atropos_n4_shrink_factor \
            -c ${class_mapping.size()} \
            -p ${basename}_%0${width}d.nii.gz \
            -q [\$(echo "\${spacing[0]}*$params.atropos_n4_bspline_fitting" | bc)] \
            -o atropos_ \
            -w $params.atropos_prior_weight

        fi

        mv atropos_Segmentation.nii.gz ${sid}_segmentation.nii.gz

        for ix in \${!tissues[@]}
        do

        mv atropos_SegmentationPosteriors\${labels[ix]}.nii.gz ${sid}_\${tissues[ix]}_pvf.nii.gz

        done

        """
}
