#!/usr/bin/env nextflow

// Diffusion Input options

// all-in-one option
params.input = false

// option for both
params.nowarp = false

params.help = false
random_generator_list = params.random_nb_generator.split(',').collect{it as int}

if(params.help) {
    usage = file("$baseDir/USAGE")

    engine = new groovy.text.SimpleTemplateEngine()

    bindings = ["flip_to_lps": "$params.flip_to_lps",
                "rois_opening": "$params.rois_opening",
                "rois_closing": "$params.rois_closing",
                "rois_smoothing": "$params.rois_smoothing",
                "rois_params": "$params.rois_params",
                "rois_seeding": "$params.rois_seeding",
                "rois_seed_params": "$params.rois_seed_params",
                "atlas": "$params.atlas",
                "flow_masked_indices": "$params.flow_masked_indices",
                "seed_masked_indices": "$params.seed_masked_indices",
                "intersections_masked_indices": "$params.intersections_masked_indices",
                "unused_connectivity_labels": "$params.unused_connectivity_labels",
                "rois_indices": "$params.rois_indices",
                "surf_smooth_nb_step": "$params.surf_smooth_nb_step",
                "surf_smooth_step_size": "$params.surf_smooth_step_size",
                "surf_flow_nb_step": "$params.surf_flow_nb_step",
                "surf_flow_step_size": "$params.surf_flow_step_size",
                "subsample_flow": "$params.subsample_flow",
                "gaussian_threshold": "$params.gaussian_threshold",
                "angle_threshold": "$params.angle_threshold",
                "random_nb_generator": "$params.random_nb_generator",
                "nb_dynamic_seeding_iter": "$params.nb_dynamic_seeding_iter",
                "seeds_weighted_per_area": "$params.seeds_weighted_per_area",
                "nb_seeds_per_random_nb": "$params.nb_seeds_per_random_nb",
                "use_seed_direction": "$params.use_seed_direction",
                "use_only_first_cut": "$params.use_only_first_cut",
                "tractography_algo": "$params.tractography_algo",
                "tractography_step": "$params.tractography_step",
                "tractography_theta": "$params.tractography_theta",
                "tractography_sfthres": "$params.tractography_sfthres",
                "pft_sfthres_init": "$params.pft_sfthres_init",
                "pft_particles": "$params.pft_particles",
                "pft_back": "$params.pft_back",
                "pft_front": "$params.pft_front",
                "compression_rate": "$params.compression_rate",
                "minimum_length": "$params.minimum_length",
                "maximum_length": "$params.maximum_length",
                "output_dir": "$params.output_dir",
                "processes": "$params.processes"]
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info ""
log.info "SET pipeline"
log.info "==============================================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if(!params.is_freesurfer && !params.is_civet){
    if (params.is_vtk){
        log.error \
"""
You need to chose another profile with "-profile vtk" to use specific atlas (i.e. "-profile vtk,freesurfer_basic")
  vtk,freesurfer_basic  :           Desikan-Killiany FS surfaces atlas with brainstem and grouped subcortical structures
  vtk,freesurfer_proper  :          Desikan-Killiany FS surfaces atlas with each subcortical structures (separate) and brainstem
  vtk,freesurfer_a2009s_basic  :    Destrieux a2009s FS surfaces atlas with brainstem and subcortical structures (together)
  vtk,freesurfer_a2009s_proper  :   Destrieux a2009s FS surfaces atlas with each subcortical structures (separate) and brainstem
"""
    }
    else{
        log.error \
"""
You need to chose a profile with "-profile" for surfaces loading
  freesurfer_basic  :           Desikan-Killiany FS surfaces atlas with brainstem and grouped subcortical structures
  freesurfer_proper  :          Desikan-Killiany FS surfaces atlas with each subcortical structures (separate) and brainstem
  freesurfer_a2009s_basic  :    Destrieux a2009s FS surfaces atlas with brainstem and subcortical structures (together)
  freesurfer_a2009s_proper  :   Destrieux a2009s FS surfaces atlas with each subcortical structures (separate) and brainstem

  civet2_dkt                :   CIVET 2.0 DKT surfaces with each subcortical structures (separate) and brainstem
  civet2_aal                :   CIVET 2.0 AAL surfaces with each subcortical structures (separate) and brainstem
"""
    }
}

process README {
    cpus 1
    tag = {"global"}
    publishDir = {"${params.output_dir}/${task.process}"}

    output:
    file "readme.txt"

    script:
    String list_options = new String();
    for (String item : params) {
        list_options += item + "\n"
    }
    """
    echo "TractoFlow pipeline\n" >> readme.txt
    echo "Start time: ${workflow.start}\n" >> readme.txt
    echo "[Command-line]\n${workflow.commandLine}\n" >> readme.txt
    echo "[Git Info]\n" >> readme.txt
    echo "${workflow.repository} - ${workflow.revision} [${workflow.commitId}]\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "${list_options}" >> readme.txt
    """
}

if (params.input){

    log.info "Input folder: ${params.input}"
    input = file(params.input)

    map_for_rois_seed = Channel
        .fromFilePairs("${input}/**/*fa.nii.gz",
                       size: 1, maxDepth:3, flat: true) {it.parent.name}
        .ifEmpty { exit 1, "Cannot find ${input}/**/*fa.nii.gz"}

    if (params.surfaceflow_only){
        tractogram_for_intersections = Channel
            .fromFilePairs("${input}/**/*.trk",
                           size: 1,
                           maxDepth:3,
                           flat: true) {it.parent.name}
    }
    else{
        fodf_and_map_for_pft = Channel
            .fromFilePairs("${input}/**/*{fodf.nii.gz,map_exclude.nii.gz,map_include.nii.gz}",
                           size: 3, maxDepth:3, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find ${input}/**/*{fodf.nii.gz,map_exclude.nii.gz,map_include.nii.gz}"}
    }

    if (!params.nowarp){
        ants_transfo_to_convert = Channel
            .fromFilePairs("${input}/**/*{output0GenericAffine.mat,output1InverseWarp.nii.gz}",
                           size: 2, maxDepth:3, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find ${input}/**/*{output0GenericAffine.mat,output1InverseWarp.nii.gz}"}
    }

    nb_sub_fodf = file("${input}/**/*fodf.nii.gz").size()
    println("Number of subject is " + nb_sub_fodf.toString())
}

if (params.is_freesurfer && params.is_civet)
{
    log.error " cannot use civet together with freesurfer profile )"
}
else if (params.is_freesurfer) {
    log.info " Profile: ${params.atlas}"
    input = file(params.input)

    if (params.atlas=="freesurfer_standard"){
        in_surfaces_label = Channel
            .fromFilePairs("${input}/**/{lh.aparc.annot,rh.aparc.annot}",
                           size: 2, maxDepth:3, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find freesurfer data: ${surfaces}/**/{lh.aparc.annot,rh.aparc.annot}" }
    }
    else if (params.atlas=="freesurfer_a2009s"){
        in_surfaces_label = Channel
            .fromFilePairs("${input}/**/{lh.aparc.a2009s.annot,rh.aparc.a2009s.annot}",
                           size: 2, maxDepth:3, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find freesurfer data ${surfaces}/**/{lh.aparc.a2009s.annot,rh.aparc.a2009s.annot}" }
    }

    in_surfaces_wmparc = Channel
        .fromFilePairs("${input}/**/wmparc*",
                       size: 1, maxDepth:3, flat: true) {it.parent.name}
        .ifEmpty { exit 1, "Cannot find freesurfer data ${surfaces}/**/wmparc*" }


    if (params.is_vtk) {
        in_surfaces_mesh = Channel
            .fromFilePairs("${input}/**/{lh*pial.vtk,lh*white.vtk,rh*pial.vtk,rh*white.vtk}",
                           size: 4, maxDepth:3, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find freesurfer data ${surfaces}/**/{lh*pial.vtk,lh*white.vtk,rh*pial.vtk,rh*white.vtk}" }


        nb_sub_surf = file("${input}/**/*lh*white.vtk").size()
        println("Number of vtk freesurfer is " + nb_sub_surf.toString())

        in_surfaces_label
            .join(in_surfaces_wmparc)
            .join(in_surfaces_mesh)
            .set{in_surfaces}

        (annots_for_surfaces_masks, annots_for_surfaces_labels, label_vol_to_convert, surfaces_for_surfaces_masks, surfaces_for_surfaces_labels, surfaces_for_lps) = in_surfaces
          .map{sid, lh_annot, rh_annot, wmparc, lh_pial, lh_white, rh_pial, rh_white ->
              [tuple(sid, lh_annot, rh_annot),
              tuple(sid, lh_annot, rh_annot),
              tuple(sid, wmparc),
              tuple(sid, lh_white, rh_white),
              tuple(sid, lh_white, rh_white),
              tuple(sid, lh_pial, rh_pial, lh_white, rh_white)]}
          .separate(6)

    }
    else{
        in_surfaces_mesh = Channel
            .fromFilePairs("${input}/**/{lh.pial,lh.white,rh.pial,rh.white}",
                           size: 4, maxDepth:3, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find freesurfer data ${surfaces}/**/{lh.pial,lh.white,rh.pial,rh.white}" }

        in_surfaces_label
            .join(in_surfaces_wmparc)
            .join(in_surfaces_mesh)
            .set{in_surfaces}

        nb_sub_surf = file("${input}/**/lh.white").size()
        println("Number of freesurfer is " + nb_sub_surf.toString())

        (annots_for_surfaces_masks, annots_for_surfaces_labels, label_vol_to_convert, freesurfer_surfaces_to_convert) = in_surfaces
          .map{sid, lh_annot, rh_annot, wmparc, lh_pial, lh_white, rh_pial, rh_white ->
              [tuple(sid, lh_annot, rh_annot),
              tuple(sid, lh_annot, rh_annot),
              tuple(sid, wmparc),
              tuple(sid, lh_pial, rh_pial, lh_white, rh_white)]}
          .separate(4)

        process A__Convert_Freesurfer_Surface {
            cpus 1

            input:
            set sid, file(lh_pial),  file(rh_pial),  file(lh_white),  file(rh_white)\
              from freesurfer_surfaces_to_convert

            output:
            set sid, "${sid}__lh_white.vtk", "${sid}__rh_white.vtk"\
              into surfaces_for_surfaces_masks, surfaces_for_surfaces_labels
            set sid, "${sid}__lh_pial.vtk", "${sid}__rh_pial.vtk",\
              "${sid}__lh_white.vtk", "${sid}__rh_white.vtk"\
              into surfaces_for_lps

            script:
            """
            mris_convert --to-scanner ${lh_pial} lh.pial.vtk
            mris_convert --to-scanner ${rh_pial} rh.pial.vtk
            mris_convert --to-scanner ${lh_white} lh.white.vtk
            mris_convert --to-scanner ${rh_white} rh.white.vtk
            mv lh.pial.vtk ${sid}__lh_pial.vtk
            mv rh.pial.vtk ${sid}__rh_pial.vtk
            mv lh.white.vtk ${sid}__lh_white.vtk
            mv rh.white.vtk ${sid}__rh_white.vtk
            """
        }
    }
}
else if (params.is_civet) {
    log.info "Input Civet: ${params.input}"
    log.info "Civet Labels: ${params.civet_template}"
    log.info " Profile: ${params.atlas}"

    if (params.atlas=="civet2"){
        input = file(params.input)
        in_civet_surf = Channel
            .fromFilePairs("${input}/**/*{gray_surface_left_81920.obj,gray_surface_right_81920.obj,white_surface_left_81920.obj,white_surface_right_81920.obj}",
                           size: 4, maxDepth:4, flat: true) {it.parent.name}
            .ifEmpty { exit 1, "Cannot find civet data: ${surfaces}/**/*{gray_surface_left_81920.obj,gray_surface_right_81920.obj,white_surface_left_81920.obj,white_surface_right_81920.obj}" }
        in_civet_transfo = Channel
            .fromFilePairs("${input}/**/*t1_tal.xfm",
                            size: 1, maxDepth:5, flat: true) {it.parent.name}
             .ifEmpty { exit 1, "Cannot find civet data: ${surfaces}/**/*t1_tal.xfm" }
        in_civet_animal = Channel
            .fromFilePairs("${input}/**/*animal_labels.mnc",
                            size: 1, maxDepth:5, flat: true) {it.parent.name}
             .ifEmpty { exit 1, "Cannot find civet data: ${surfaces}/**/*animal_labels.mnc" }
    }

    process A__Civet_Template {
        cpus 1
        tag = {"global"}
        publishDir = {"${params.output_dir}/${task.process}"}

        output:
        set "template_left.txt", "template_right.txt"\
            into in_civet_template

        script:
        """
        cp ${file(params.civet_template)}/*left.txt template_left.txt
        cp ${file(params.civet_template)}/*right.txt template_right.txt
        """
    }

    in_civet_surf
        .join(in_civet_transfo)
        .join(in_civet_animal)
        .combine(in_civet_template)
        .set{in_civet}

    nb_sub_surf = file("${input}/**/*gray_surface_left_81920.obj").size()
    println("Number of civet is " + nb_sub_surf.toString())

    (annots_for_surfaces_masks, annots_for_surfaces_labels, xfm_transfo_to_convert, animal_to_convert, civet_surfaces_to_convert) = in_civet
      .map{sid, lh_pial, rh_pial, lh_white, rh_white, xfm_transfo, animal_labels, lh_annot, rh_annot  ->
          [tuple(sid, lh_annot, rh_annot),
          tuple(sid, lh_annot, rh_annot),
          tuple(sid, xfm_transfo),
          tuple(sid, animal_labels, xfm_transfo),
          tuple(sid, lh_pial, rh_pial, lh_white, rh_white, xfm_transfo)]}
      .separate(5)

    process A__Convert_CIVET_Surface {
        cpus 1

        input:
        set sid, file(lh_pial), file(rh_pial), file(lh_white), file(rh_white), file(xfm_transfo)\
            from civet_surfaces_to_convert

        output:
        set sid, "${sid}__lh_white.vtk", "${sid}__rh_white.vtk"\
            into surfaces_for_surfaces_masks, surfaces_for_surfaces_labels
        set sid, "${sid}__lh_pial.vtk", "${sid}__rh_pial.vtk",\
            "${sid}__lh_white.vtk", "${sid}__rh_white.vtk"\
            into surfaces_for_lps

        script:
        """
        xfminvert ${xfm_transfo}  ${sid}__t1_tal_inv.xfm
        sed '''/Linear_Transform =/,/;/!d ; /Linear/d ; s/^ //; s/;/\\n0 0 0 1/g'''\
            ${xfm_transfo} > ${sid}__to_t1_transfo.txt
        sed '''/Linear_Transform =/,/;/!d ; /Linear/d ; s/^ //; s/;/\\n0 0 0 1/g'''\
            ${sid}__t1_tal_inv.xfm > ${sid}__to_t1_inv_transfo.txt

        transform_objects ${lh_pial} ${sid}__t1_tal_inv.xfm lh_pial_t1.mni.obj
        transform_objects ${rh_pial} ${sid}__t1_tal_inv.xfm rh_pial_t1.mni.obj
        transform_objects ${lh_white} ${sid}__t1_tal_inv.xfm lh_white_t1.mni.obj
        transform_objects ${rh_white} ${sid}__t1_tal_inv.xfm rh_white_t1.mni.obj

        scil_convert_surface.py lh_pial_t1.mni.obj ${sid}__lh_pial.vtk
        scil_convert_surface.py rh_pial_t1.mni.obj ${sid}__rh_pial.vtk
        scil_convert_surface.py lh_white_t1.mni.obj ${sid}__lh_white.vtk
        scil_convert_surface.py rh_white_t1.mni.obj ${sid}__rh_white.vtk
        """
    }

    process A__Convert_Animal {
        cpus 1

        input:
        set sid, file(animal_labels), file(xfm_transfo)\
            from animal_to_convert

        output:
        set sid, "${sid}__animal_labels_native.nii"\
            into label_vol_to_convert

        script:
        """
        mincresample -nearest_neighbour -tfm_input_sampling\
            -invert_transformation -transformation ${xfm_transfo}\
            ${animal_labels} ${sid}__animal_labels_native.mnc

        mnc2nii ${sid}__animal_labels_native.mnc ${sid}__animal_labels_native.nii
        """
    }
}
else{
    log.error "Use a profile (freesurfer: -profile freesurfer_proper , civet: -profile civet2_dkt)"
}

// setup variable
nb_subject = [nb_sub_fodf, nb_sub_surf].min()
println("Number of subject (with surface and fodf) is " + nb_subject.toString())

process A__Convert_Label_Volume {
    cpus 1

    input:
    set sid, file(label_vol)\
        from label_vol_to_convert

    output:
    set sid, "${sid}__labels.nii.gz"\
        into label_vol_for_rois

    script:
    """
    mri_convert ${label_vol} ${sid}__labels.nii.gz
    """
}

process A__Surface_to_LPS{
    cpus 1

    input:
    set sid, file(lh_pial), file(rh_pial), file(lh_white), file(rh_white)\
        from surfaces_for_lps

    output:
    set sid, "${sid}__lh_pial_lps.vtk", "${sid}__rh_pial_lps.vtk",\
        "${sid}__lh_white_lps.vtk", "${sid}__rh_white_lps.vtk"\
        into surfaces_to_concatenate

    script:
    """
    scil_flip_surface.py ${lh_pial} ${sid}__lh_pial_lps.vtk ${params.flip_to_lps}
    scil_flip_surface.py ${rh_pial} ${sid}__rh_pial_lps.vtk ${params.flip_to_lps}
    scil_flip_surface.py ${lh_white} ${sid}__lh_white_lps.vtk ${params.flip_to_lps}
    scil_flip_surface.py ${rh_white} ${sid}__rh_white_lps.vtk ${params.flip_to_lps}
    """
}

annots_for_surfaces_masks
    .join(surfaces_for_surfaces_masks)
    .set{data_for_surfaces_masks}

process B__Surface_Mask {
    cpus 1

    input:
    set sid, file(lh_annot), file(rh_annot), file(lh_surf), file(rh_surf)\
        from data_for_surfaces_masks

    output:
    set sid,  "${sid}__lh_flow_mask.npy", "${sid}__rh_flow_mask.npy",\
        "${sid}__lh_seed_mask.npy", "${sid}__rh_seed_mask.npy",\
        "${sid}__lh_intersections_mask.npy", "${sid}__rh_intersections_mask.npy",\
        "${sid}__lh_zero_mask.npy", "${sid}__rh_zero_mask.npy"\
        into masks_for_concatenate

    script:
    label_tag="--annot"
    if (params.is_civet) {
        label_tag="--vts_label"
    }
    """
    scil_surface.py ${lh_surf} ${label_tag} ${lh_annot}\
        -i ${params.flow_masked_indices} --inverse_mask\
        --save_vts_mask ${sid}__lh_flow_mask.npy
    scil_surface.py ${rh_surf} ${label_tag} ${rh_annot}\
        -i ${params.flow_masked_indices} --inverse_mask\
        --save_vts_mask ${sid}__rh_flow_mask.npy

    scil_surface.py ${lh_surf} ${label_tag} ${lh_annot}\
        -i ${params.seed_masked_indices} --inverse_mask\
        --save_vts_mask ${sid}__lh_seed_mask.npy
    scil_surface.py ${rh_surf} ${label_tag} ${rh_annot}\
        -i ${params.seed_masked_indices} --inverse_mask\
        --save_vts_mask ${sid}__rh_seed_mask.npy

    scil_surface.py ${lh_surf} ${label_tag} ${lh_annot}\
        -i ${params.intersections_masked_indices} --inverse_mask\
        --save_vts_mask ${sid}__lh_intersections_mask.npy
    scil_surface.py ${rh_surf} ${label_tag} ${rh_annot}\
        -i ${params.intersections_masked_indices} --inverse_mask\
        --save_vts_mask ${sid}__rh_intersections_mask.npy

    scil_surface.py ${lh_surf} --vts_val 0.0 --save_vts_mask ${sid}__lh_zero_mask.npy
    scil_surface.py ${rh_surf} --vts_val 0.0 --save_vts_mask ${sid}__rh_zero_mask.npy
    """
}

surfaces_for_surfaces_labels
    .join(annots_for_surfaces_labels)
    .set{data_for_surfaces_labels}

process B__Surface_Label {
    cpus 1

    input:
    set sid, file(lh_surf), file(rh_surf), file(lh_annot), file(rh_annot)\
        from data_for_surfaces_labels

    output:
    set sid, "${sid}__lh_labels.npy", "${sid}__rh_labels.npy", "${sid}__lh_zero_mask.npy", "${sid}__rh_zero_mask.npy"\
        into labels_for_cocatenate

    script:
    label_tag=" --annot "
    if (params.is_civet) {
        label_tag=" --vts_label "
    }
    """
    scil_surface.py ${lh_surf} ${label_tag} ${lh_annot} --save_vts_label ${sid}__lh_labels.npy
    scil_surface.py ${rh_surf} ${label_tag} ${rh_annot} --save_vts_label ${sid}__rh_labels.npy
    scil_surface.py ${lh_surf} --vts_val 0.0 --save_vts_mask ${sid}__lh_zero_mask.npy
    scil_surface.py ${rh_surf} --vts_val 0.0 --save_vts_mask ${sid}__rh_zero_mask.npy
    """
}

// Generate ROIs in VTK
process B__Generate_ROI {
    cpus 1

    input:
    set sid, file(vol)\
        from label_vol_for_rois

    output:
    set sid, "${sid}__roi*.vtk"\
        into rois_to_concatenate, rois_for_masks

    script:
    command_lines=""
    params.rois_indices.eachWithIndex{ roiString, index ->
        command_lines +=\
            """
            scil_surface_from_volume.py $vol\
                ${sid}__roi${index.toString().padLeft(4, "0")}.vtk\
                --index ${roiString}\
                --closing ${params.rois_closing}\
                --opening ${params.rois_opening}\
                --smooth ${params.rois_smoothing}\
                --vox2vtk ${params.rois_params}
            """
    }
    """
    $command_lines
    """
}

map_for_rois_seed
    .join(rois_for_masks)
    .set{map_and_rois_for_masks}

process B__ROI_Mask {
    cpus 1

    input:
    set sid, file(map), file(rois)\
        from map_and_rois_for_masks

    output:
    set sid, "${sid}__roi*_flow_mask.npy", "${sid}__roi*_seed_mask.npy", "${sid}__roi*_intersections_mask.npy"\
        into rois_mask_for_concatenate

    set sid, "${sid}__roi*_intersections_mask.npy"\
        into rois_labels_for_concatenate

    when:
    params.rois_seeding

    script:
    command_lines=""
    rois.each{
        if (params.rois_seeding)
        {
            command_lines += \
                """
                scil_surface.py ${it} --vts_val 0.0 --save_vts_mask ${it.getSimpleName()}_flow_mask.npy \n
                scil_surface.py ${it} --vts_val 1.0 --save_vts_mask ${it.getSimpleName()}_intersections_mask.npy \n
                scil_surface_map_from_volume.py ${it} ${map} ${it.getSimpleName()}_seed_mask.npy ${params.rois_seed_params} \n
                """
        }
        else{
            command_lines += \
                """
                scil_surface.py ${it} --vts_val 0.0 --save_vts_mask ${it.getSimpleName()}_seed_mask.npy \n
                scil_surface.py ${it} --vts_val 0.0 --save_vts_mask ${it.getSimpleName()}_flow_mask.npy \n
                scil_surface.py ${it} --vts_val 1.0 --save_vts_mask ${it.getSimpleName()}_intersections_mask.npy \n
                """
        }
    }
    """
    ${command_lines}
    """
}


surfaces_to_concatenate
    .join(rois_to_concatenate)
    .set{surfaces_and_rois_to_concatenate}

process B__Concatenate_Surface {
    cpus 1

    input:
    set sid, file(lh_pial), file(rh_pial), file(lh_white), file(rh_white), file(rois)\
        from surfaces_and_rois_to_concatenate

    output:
    set sid, "${sid}__surfaces.vtk"\
        into surfaces_to_warp
    set sid, "${sid}__surfaces_type.npy"\
        into surface_type_for_set_nf

    file("${sid}__surfaces_id.npy")

    script:
    rois_list=""
    rois.each{
        rois_list +="${it} "
    }
    """
    scil_concatenate_surfaces.py ${lh_white} ${rh_white}\
        --outer_surfaces ${lh_pial} ${rh_pial}\
        --inner_surfaces ${rois_list}\
        --out_surface_id ${sid}__surfaces_id.npy\
        --out_surface_type_map ${sid}__surfaces_type.npy\
        --out_concatenated_surface ${sid}__surfaces.vtk
    """
}

masks_for_concatenate
    .join(rois_mask_for_concatenate)
    .set{all_masks_for_concatenate}

process B__Concatenate_Mask {
    cpus 1

    input:
    set sid, file(lh_flow_mask), file(rh_flow_mask), file(lh_seed_mask), file(rh_seed_mask),\
        file(lh_intersections_mask), file(rh_intersections_mask), file(lh_zero_mask), file(rh_zero_mask),\
        file(rois_flow_m), file(rois_seed_m), file(rois_intersections_m)\
            from all_masks_for_concatenate

    output:
    set sid, "${sid}__flow_mask.npy"\
        into flow_mask_for_set_nf
    set sid, "${sid}__seed_mask.npy"\
        into seed_mask_for_set_nf
    set sid, "${sid}__intersections_mask.npy"\
        into intersections_mask_for_set_nf

    script:
    rois_flow_masks=""
    rois_seed_masks=""
    rois_intersections_masks=""

    rois_flow_m.each{rois_flow_masks +="${it} "}
    rois_seed_m.each{rois_seed_masks +="${it} "}
    rois_intersections_m.each{rois_intersections_masks +="${it} "}

    """
    scil_concatenate_surfaces_map.py ${lh_flow_mask} ${rh_flow_mask}\
        --outer_surfaces_map ${lh_zero_mask} ${rh_zero_mask}\
        --inner_surfaces_map ${rois_flow_m}\
        --out_map ${sid}__flow_mask.npy
    scil_concatenate_surfaces_map.py ${lh_seed_mask} ${rh_seed_mask}\
        --outer_surfaces_map ${lh_zero_mask} ${rh_zero_mask}\
        --inner_surfaces_map ${rois_seed_m}\
        --out_map ${sid}__seed_mask.npy
    scil_concatenate_surfaces_map.py ${lh_intersections_mask} ${rh_intersections_mask}\
        --outer_surfaces_map ${lh_intersections_mask} ${rh_intersections_mask}\
        --inner_surfaces_map ${rois_intersections_m}\
        --out_map ${sid}__intersections_mask.npy
    """
}

labels_for_cocatenate
    .join(rois_labels_for_concatenate)
    .set{all_labels_for_concatenate}

process B__Concatenate_Label {
    cpus 1

    input:
    set sid, file(lh_labels), file(rh_labels), file(lh_zero), file(rh_zero), file(rois_labels)\
        from all_labels_for_concatenate

    output:
    set sid, "${sid}__unique_id.npy"\
        into labels_for_set_nf

    file("${sid}__unique_id.txt")

    script:
    rois_labels_s=""
    rois_labels.each{rois_labels_s +="${it} "}
    unused_labels_tag=""
    if (params.unused_connectivity_labels && params.unused_connectivity_labels?.trim()) {
        unused_labels_tag=" --indices_to_remove ${params.unused_connectivity_labels}"
    }
    """
    scil_concatenate_surfaces_map.py ${lh_labels} ${rh_labels} \
        --outer_surfaces_map ${lh_zero} ${rh_zero}\
        --inner_surfaces_map ${rois_labels_s}\
        --out_map ${sid}__unique_id.npy\
        --unique_id\
        --out_id_map ${sid}__unique_id.txt\
        ${unused_labels_tag}
    """
}

// Transform surfaces and ROIs
if (params.nowarp){
    surfaces_to_warp
        .into{surfaces_for_seed;surfaces_for_smooth;surfaces_for_density}
}
else{
    // if "tractoflow" or "antswarp"
    process A__Convert_ANTs_Transformation {
        cpus 1

        input:
        set sid, file(affine_transfo), file(warp_transfo)\
            from ants_transfo_to_convert

        output:
        set sid, "${sid}__vtk_transfo.txt", file(warp_transfo)\
            into ants_transfo_for_surfaces

        script:
        """
        ConvertTransformFile 3 ${affine_transfo}\
            ${sid}__vtk_transfo.txt --hm
        """
    }

    ants_transfo_for_surfaces
        .join(surfaces_to_warp)
        .set{surfaces_with_transform}

    //Apply Transform to surfaces (from t1 to dwi space)
    process C__Register_Surface {
        cpus 1

        input:
        set sid, file(affine_transfo), file(warp_transfo), file(surfaces)\
            from surfaces_with_transform

        output:
        set sid, "${sid}__surfaces_b0.vtk"\
            into surfaces_for_seed, surfaces_for_smooth, surfaces_for_density

        script:
        """
        scil_transform_surface.py ${surfaces} ${affine_transfo}\
            ${sid}__surfaces_b0.vtk\
            --ants_warp ${warp_transfo}
        """
    }
}

surfaces_for_smooth
    .join(flow_mask_for_set_nf)
    .set{data_for_surface_smooth}

process D__Surface_Flow {
    cpus 1

    input:
    set sid, file(surf), file(mask)\
        from data_for_surface_smooth

    output:
    set sid, "${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.vtk"\
        into surface_flow_surfaces_for_pft, surfaces_for_connectivity
    set sid, "${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.hdf5"\
        into surface_flow_lines_for_combine

    file("${sid}__smoothed.vtk")

    script:
    if ((params.surf_flow_nb_step as Integer) > 1 )
        """
        scil_smooth_surface.py ${surf} ${sid}__smoothed.vtk\
            --vts_mask ${mask}\
            --nb_steps ${params.surf_smooth_nb_step}\
            --step_size ${params.surf_smooth_step_size}\

        scil_surface_flow.py ${sid}__smoothed.vtk\
            ${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.vtk\
            --vts_mask ${mask}\
            --nb_step ${params.surf_flow_nb_step}\
            --step_size ${params.surf_flow_step_size}\
            --subsample_flow ${params.subsample_flow}\
            --gaussian_threshold ${params.gaussian_threshold}\
            --angle_threshold ${params.angle_threshold}\
            --out_flow ${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.hdf5
        """
    else
        """
        scil_smooth_surface.py ${surf} ${sid}__smoothed.vtk\
            --vts_mask ${mask}\
            --nb_steps ${params.surf_smooth_nb_step}\
            --step_size ${params.surf_smooth_step_size}\

        cp ${sid}__smoothed.vtk ${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.vtk
        touch ${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.hdf5
        """
}


if(params.surfaceflow_only){
    surface_flow_surfaces_for_pft
    .combine(intersections_mask_for_set_nf, by : 0)
    .combine(surface_type_for_set_nf, by : 0)
    .combine(surface_flow_lines_for_combine, by : 0)
        .set{surface_and_map_for_intersections}

    tractogram_for_intersections
        .combine(surface_and_map_for_intersections, by : 0)
        .set{data_for_intersections}


    process EF__Surface_Flow_onto_Tractography {
        cpus 1

        input:
        set sid, file(tractogram), file(surf), file(s_mask), file(s_type), file(flow)\
                from data_for_intersections
        // set sid, file(tractogram), file(surf), file(flow), file(s_mask), file(s_type)\
        //     from data_for_intersections
        // set sid, file(seed_map), file(sum_density), random_id, loop_id,\
        //     file(surf), file(fodf), file(map_exclude), file(map_include),\
        //     file(s_mask), file(s_type), file(flow)\
        //         from data_for_set

        output:
        set sid, "${tractogram.getSimpleName()}_filtered.npz"\
            into intersections_for_concatenate
        set sid, "${tractogram.getSimpleName()}_filtered.fib"\
            into streamlines_for_concatenate

        script:
        surf_flow_command="cp ${tractogram.getSimpleName()}_cut.fib ${tractogram.getSimpleName()}_sf.fib"
        if ((params.surf_flow_nb_step as Integer) > 1 )
            surf_flow_command=\
                """
                scil_surface_combine_flow.py ${surf} \
                    ${flow} \
                    ${tractogram.getSimpleName()}_cut.npz \
                    ${tractogram.getSimpleName()}_cut.fib \
                    ${tractogram.getSimpleName()}_sf.fib \
                    --compression_rate ${params.compression_rate}
                """
        """
        scil_convert_tractogram.py ${tractogram} ${tractogram.getSimpleName()}.fib

        scil_surface_tractogram_intersections.py ${surf} \
            ${tractogram.getSimpleName()}.fib \
            ${s_type} ${s_mask} \
            --output_intersections ${tractogram.getSimpleName()}_cut.npz \
            --output_tractogram ${tractogram.getSimpleName()}_cut.fib

        $surf_flow_command

        scil_surface_filtering.py  ${surf}\
            ${tractogram.getSimpleName()}_cut.npz\
            ${tractogram.getSimpleName()}_sf.fib \
            ${tractogram.getSimpleName()}_filtered.fib\
            --out_intersection ${tractogram.getSimpleName()}_filtered.npz\
            --min_length ${params.minimum_length}\
            --max_length ${params.maximum_length}

        """
    }

    intersections_for_concatenate
        .groupTuple()
        .set{intersections_grouped_for_concatenate}

}
else{

    surfaces_for_seed
        .join(seed_mask_for_set_nf)
        .set{data_rand_for_surface_seed}

    process E__Surface_Seeding_Map {
        cpus 1

        input:
        set sid, file(surf), file(mask)\
            from data_rand_for_surface_seed

        output:
        set sid, "${sid}__seeding_map_0.npy", "${sid}__zeros_tri_map.npy"\
            into surfaces_seeding_map_for_set

        script:
        area_seeding_tag=""
        if (params.seeds_weighted_per_area) {
            area_seeding_tag=" --triangle_area_weighting "
        }
        """
        scil_surface_seed_map.py ${surf}\
            ${sid}__seeding_map_0.npy\
            --vts_mask ${mask}\
            ${area_seeding_tag}

        scil_surface_seed_map.py ${surf} ${sid}__zeros_tri_map.npy\
            --zeros_map
        """
    }

    // Where the magic happen
    nb_tracking_per_sub = params.nb_dynamic_seeding_iter * (random_generator_list).size()
    total_tracking = nb_subject*(nb_tracking_per_sub+1)

    surfaces_seeding_map_for_set
        .combine(random_generator_list)
        .combine([0])
        .set{setup_loop_ch}

    feedback_ch = Channel.create()
    intersections_for_concatenate = Channel.create()
    streamlines_for_concatenate = Channel.create()
    set_input_ch = setup_loop_ch.mix(feedback_ch).take(total_tracking)

    set_input_ch
        .combine(surface_flow_surfaces_for_pft, by : 0)
        .combine(fodf_and_map_for_pft, by : 0)
        .combine(intersections_mask_for_set_nf, by : 0)
        .combine(surface_type_for_set_nf, by : 0)
        .combine(surface_flow_lines_for_combine, by : 0)
        .set{data_for_set}

    process F__Surface_Enhanced_Tractography {
        cpus 1

        input:
        set sid, file(seed_map), file(sum_density), random_id, loop_id,\
            file(surf), file(fodf), file(map_exclude), file(map_include),\
            file(s_mask), file(s_type), file(flow)\
                from data_for_set

        output:
        set sid, "${sid}__intersections_${rand_loop_id}_filtered.npz"\
            into intersections_for_concatenate
        set sid, "${sid}__set_${rand_loop_id}_filtered.fib"\
            into streamlines_for_concatenate
        set sid, file(seed_map), "${sid}__sum_density_${rand_loop_id}.npy", random_id, next_id\
            into feedback_ch

        file "${sid}__seeding_map_${rand_loop_id}.npy"
        file "${sid}__seeds_${rand_loop_id}.npz"
        file "${sid}__set_density_${rand_loop_id}.npy"

        when:
        loop_id < params.nb_dynamic_seeding_iter

        script:
        next_id = loop_id + 1
        rand_loop_id = random_id.toString().padLeft(4, "0") + "_i" + loop_id.toString().padLeft(4, "0")

        seed_direction_tag=""
        seed_direction_tag_inter=""
        if (params.use_seed_direction) {
            seed_direction_tag="--set_dir"
            seed_direction_tag_inter="--surface_seeds ${sid}__seeds_${rand_loop_id}.npz"
        }

        first_cut_tag=""
        if (params.use_only_first_cut) {
            first_cut_tag="--only_first_cut"
        }

        flow_line="cp ${sid}__cut_${rand_loop_id}.fib ${sid}__set_${rand_loop_id}.fib"
        if ((params.surf_flow_nb_step as Integer) > 1 ) {
            flow_line=""" scil_surface_combine_flow.py ${surf} ${flow}\
                        ${sid}__intersections_${rand_loop_id}.npz\
                        ${sid}__cut_${rand_loop_id}.fib\
                        ${sid}__set_${rand_loop_id}.fib\
                        --compression_rate ${params.compression_rate} """
        }

        """
        scil_surface_seed_map.py ${surf}\
            ${sid}__seeding_map_${rand_loop_id}.npy\
            --triangle_weight ${seed_map}\
            --previous_density ${sum_density}

        scil_surface_seeds_from_map.py ${surf} ${sid}__seeding_map_${rand_loop_id}.npy\
            ${params.nb_seeds_per_random_nb}\
            ${sid}__seeds_${rand_loop_id}.npz\
            --random_number_generator ${random_id}

        scil_surface_pft_dipy.py ${fodf} ${map_include} ${map_exclude} ${surf}\
            ${sid}__seeds_${rand_loop_id}.npz\
            ${sid}__set_${rand_loop_id}.trk\
            --algo ${params.tractography_algo}\
            --step ${params.tractography_step}\
            --theta ${params.tractography_theta}\
            --sfthres ${params.tractography_sfthres}\
            --max_length ${params.maximum_length}\
            --random_seed ${loop_id}\
            --compress ${params.compression_rate}\
            --particles ${params.pft_particles}\
            --back ${params.pft_back}\
            --forward ${params.pft_front}\
            ${seed_direction_tag}

        scil_convert_tractogram.py ${sid}__set_${rand_loop_id}.trk \
            ${sid}__set_${rand_loop_id}.fib

        scil_surface_tractogram_intersections.py ${surf}\
            ${sid}__set_${rand_loop_id}.fib\
            ${s_type} ${s_mask} ${seed_direction_tag_inter} ${first_cut_tag}\
            --output_intersections ${sid}__intersections_${rand_loop_id}.npz\
            --output_tractogram ${sid}__cut_${rand_loop_id}.fib\

        $flow_line

        scil_surface_filtering.py ${surf}\
            ${sid}__intersections_${rand_loop_id}.npz\
            ${sid}__set_${rand_loop_id}.fib\
            ${sid}__set_${rand_loop_id}_filtered.fib\
            --out_intersections ${sid}__intersections_${rand_loop_id}_filtered.npz\
            --min_length ${params.minimum_length}\
            --max_length ${params.maximum_length}

        scil_surface_intersections_density.py ${surf} ${sid}__intersections_${rand_loop_id}_filtered.npz\
            ${sid}__set_density_${rand_loop_id}.npy

        scil_surface_seed_map.py ${surf} ${sid}__sum_density_${rand_loop_id}.npy\
            --sum_maps ${sum_density} ${sid}__set_density_${rand_loop_id}.npy
        """
    }

    intersections_for_concatenate
        .groupTuple(size: nb_tracking_per_sub)
        .set{intersections_grouped_for_concatenate}
}

process G__Concatenate_Intersection {
    cpus 1

    input:
    set sid, file(intersections)\
        from intersections_grouped_for_concatenate

    output:
    set sid, "${sid}__set_c_filtered.npz"\
        into intersections_for_connectivity, intersections_for_density

    script:
    command_lines=""
    intersections.each{
        command_lines +=" ${it} "
    }
    """
    scil_concatenate_surfaces_intersections.py ${command_lines}\
        --output_intersections ${sid}__set_c_filtered.npz
    """
}

surfaces_for_connectivity
    .join(labels_for_set_nf)
    .set{surf_and_labels_for_connectivity}

intersections_for_connectivity
    .combine(surf_and_labels_for_connectivity, by : 0)
    .set{data_for_connectivity}

process H__Compute_Connectivity_Matrix {
    cpus 1

    input:
    set sid, file(intersections), file(surf), file(labels)\
        from data_for_connectivity

    output:
    set sid, "${sid}__set_connectivity.npy"\
        into matrices_for_stats

    script:
    """
    scil_surface_intersections_to_connectivity.py ${surf} ${intersections} ${labels}\
        ${sid}__set_connectivity.npy
    """
}

intersections_for_density
    .combine(surfaces_for_density, by : 0)
    .set{data_for_density}

process H__Compute_Surface_Density {
    cpus 1

    input:
    set sid, file(intersections), file(surf)\
        from data_for_density

    output:
    set sid, "${sid}__set_density.npy"\
        into density_for_stats

    script:
    """
    scil_surface_intersections_density.py ${surf} ${intersections}\
    ${sid}__set_density.npy --normalize_l1_to 1
    """
}

// streamlines_for_concatenate
//     .groupTuple(size: nb_tracking_per_sub)
//     .set{streamlines_grouped_for_concatenate}
//
// process G__Concatenate_Tractogram {
//     cpus 1
//
//     input:
//     set sid, file(tractograms) from streamlines_grouped_for_concatenate
//
//     output:
//     file "${sid}__set_c_filtered.fib"
//
//     when:
//     params.concatenate_tractogram
//
//     script:
//     command_lines=""
//     tractograms.each{
//         command_lines +=" ${it} "
//     }
//     """
//     scil_streamlines_math.py concatenate $command_lines\
//         ${sid}__set_c_filtered.fib
//     """
// }
