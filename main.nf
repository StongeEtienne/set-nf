#!/usr/bin/env nextflow

// Surface input options
params.surfaces = false

// Tractogram input options
params.tractoflow = false

params.nowarp = false

params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    engine = new groovy.text.SimpleTemplateEngine()

    bindings = ["flip_to_lps": "$params.flip_to_lps",
                "rois_opening": "$params.rois_opening",
                "rois_closing": "$params.rois_closing",
                "rois_smoothing": "$params.rois_smoothing",
                "rois_params": "$params.rois_params",
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
log.info "SET FLOW pipeline"
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
  vtk,freesurfer_basic          : Freesurfer surfaces with brainstem and subcortical structures (together)
  vtk,freesurfer_proper         : Freesurfer surfaces with each subcortical structures (separate) and brainstem
  vtk,freesurfer_a2009s_basic   : Freesurfer a2009s surfaces with brainstem and subcortical structures (together)
  vtk,freesurfer_a2009s_proper  : Freesurfer a2009s surfaces with each subcortical structures (separate) and brainstem
"""
    }
    else{
        log.error \
"""
You need to chose a profile with "-profile" for surfaces loading
  freesurfer_basic          : Freesurfer surfaces with brainstem and subcortical structures (together)
  freesurfer_proper         : Freesurfer surfaces with each subcortical structures (separate) and brainstem
  freesurfer_a2009s_basic   : Freesurfer a2009s surfaces with brainstem and subcortical structures (together)
  freesurfer_a2009s_proper  : Freesurfer a2009s surfaces with each subcortical structures (separate) and brainstem

  civet2_dkt                : CIVET 2.0 DKT surfaces with each subcortical structures (separate) and brainstem
  civet2_aal                : CIVET 2.0 AAL surfaces with each subcortical structures (separate) and brainstem
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

if (params.tractoflow){
    log.info "Input: ${params.tractoflow}"
    tractoflow = file(params.tractoflow)
    tractogram_for_intersections = Channel
        .fromFilePairs("${tractoflow}/**/Tracking/*.trk",
                       size: 1,
                       maxDepth:3,
                       flat: true) {it.parent.parent.name}

    ants_transfo_to_convert = Channel
        .fromFilePairs("${tractoflow}/**/Register_T1/*{output0GenericAffine.mat,output1InverseWarp.nii.gz}",
                       size: 2,
                       maxDepth:3,
                       flat: true) {it.parent.parent.name}
}
else {
    log.error "Use --tractoflow /path/to/tractoflow/results"
}


// *****************************************
// START COPY FROM SET-NF (with minor modification)
// *****************************************
if (params.is_freesurfer && params.is_civet)
{
    log.error " cannot use civet together with freesurfer profile )"
}
else if (params.is_freesurfer) {
    log.info "Input Freesurfer: ${params.surfaces}"
    log.info " Profile: ${params.atlas}"
    surfaces = file(params.surfaces)

    if (params.is_vtk) {
        if (params.atlas=="freesurfer_standard"){
            in_surfaces = Channel
                .fromFilePairs("$surfaces/**/*{lh*aparc.annot,lh*pial.vtk,lh*white.vtk,rh*aparc.annot,rh*pial.vtk,rh*white.vtk,wmparc*}",
                               size: 7,
                               maxDepth:3,
                               flat: true) {it.parent.name}
        }
        else if (params.atlas=="freesurfer_a2009s"){
            in_surfaces = Channel
                .fromFilePairs("$surfaces/**/*{lh*a2009s.annot,lh*pial.vtk,lh*white.vtk,rh*a2009s.annot,rh*pial.vtk,rh*white.vtk,wmparc*}",
                               size: 7,
                               maxDepth:3,
                               flat: true) {it.parent.name}
        }
        else
        {
            log.error "Freesurfer profile should be given with vtk input ( e.g. --profile vtk, freesurfer_proper )"
        }

        (annots_for_surfaces_masks, annots_for_surfaces_labels, label_vol_to_convert, surfaces_for_surfaces_masks, surfaces_for_surfaces_labels, surfaces_for_lps) = in_surfaces
          .map{sid, lh_annot, lh_pial, lh_white, rh_annot, rh_pial, rh_white, wmparc ->
              [tuple(sid, lh_annot, rh_annot),
              tuple(sid, lh_annot, rh_annot),
              tuple(sid, wmparc),
              tuple(sid, lh_white, rh_white),
              tuple(sid, lh_white, rh_white),
              tuple(sid, lh_pial, rh_pial, lh_white, rh_white)]}
          .separate(6)
    }
    else{
        if (params.atlas=="freesurfer_standard"){
            in_surfaces = Channel
                .fromFilePairs("${surfaces}/**/*{label/lh.aparc.annot,label/rh.aparc.annot,mri/wmparc*,surf/lh.pial,surf/lh.white,surf/rh.pial,surf/rh.white}",
                               size: 7,
                               maxDepth:4,
                               flat: true) {it.parent.parent.name}
        }
        else if (params.atlas=="freesurfer_a2009s"){
            in_surfaces = Channel
                .fromFilePairs("${surfaces}/**/*{label/lh.aparc.a2009s.annot,label/rh.aparc.a2009s.annot,mri/wmparc*,surf/lh.pial,surf/lh.white,surf/rh.pial,surf/rh.white}",
                               size: 7,
                               maxDepth:4,
                               flat: true) {it.parent.parent.name}
        }

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
    log.info "Input Civet: ${params.surfaces}"
    log.info "Civet Labels: ${params.civet_template}"
    log.info " Profile: ${params.atlas}"

    if (params.atlas=="civet2"){
        surfaces = file(params.surfaces)
        in_civet_surf = Channel
            .fromFilePairs("${surfaces}/**/{surfaces/*gray_surface_left_81920.obj,surfaces/*gray_surface_right_81920.obj,surfaces/*white_surface_left_81920.obj,surfaces/*white_surface_right_81920.obj}",
                           size: 4,
                           maxDepth:4,
                           flat: true) {it.parent.parent.name}
        in_civet_transfo = Channel
            .fromFilePairs("${surfaces}/**/transforms/linear/*t1_tal.xfm",
                            size: 1,
                            maxDepth:5,
                            flat: true) {it.parent.parent.parent.name}
        in_civet_animal = Channel
            .fromFilePairs("${surfaces}/**/segment/*animal_labels.mnc",
                            size: 1,
                            maxDepth:5,
                            flat: true) {it.parent.parent.name}
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
        set sid, file(lh_pial),  file(rh_pial),  file(lh_white),  file(rh_white), file(xfm_transfo)\
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
    log.error "Use a profile (freesurfer: --profile freesurfer_proper , civet: --profile civet2_dkt)"
}

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
        into labels_for_concatenate

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

process B__ROI_Mask {
    cpus 1

    input:
    set sid, file(rois)\
        from rois_for_masks

    output:
    set sid, "${sid}__roi*_flow_mask.npy", "${sid}__roi*_seed_mask.npy", "${sid}__roi*_intersections_mask.npy"\
        into rois_mask_for_concatenate

    set sid, "${sid}__roi*_intersections_mask.npy"\
        into rois_labels_for_concatenate

    script:
    command_lines=""
    rois.each{
        command_lines +="scil_surface.py ${it} --vts_val 0.0 --save_vts_mask ${it.getSimpleName()}_seed_mask.npy \n"
        command_lines +="scil_surface.py ${it} --vts_val 0.0 --save_vts_mask ${it.getSimpleName()}_flow_mask.npy \n"
        command_lines +="scil_surface.py ${it} --vts_val 1.0 --save_vts_mask ${it.getSimpleName()}_intersections_mask.npy \n"
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
        into surface_type_for_sf

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
        into flow_mask_for_sf
    set sid, "${sid}__seed_mask.npy"\
        into seed_mask_for_sf
    set sid, "${sid}__intersections_mask.npy"\
        into intersections_mask_for_sf

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

labels_for_concatenate
    .join(rois_labels_for_concatenate)
    .set{all_labels_for_concatenate}

process B__Concatenate_Label {
    cpus 1

    input:
    set sid, file(lh_labels), file(rh_labels), file(lh_zero), file(rh_zero), file(rois_labels)\
        from all_labels_for_concatenate

    output:
    set sid, "${sid}__unique_id.npy"\
        into labels_for_sf

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
// REMOVED SEED
if (params.nowarp){
    surfaces_to_warp
        .into{surfaces_for_smooth;surfaces_for_connectivity;surfaces_for_density}
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
    // REMOVED SEED
    process C__Register_Surface {
        cpus 1

        input:
        set sid, file(affine_transfo), file(warp_transfo), file(surfaces)\
            from surfaces_with_transform

        output:
        set sid, "${sid}__surfaces_b0.vtk"\
            into surfaces_for_smooth, surfaces_for_connectivity, surfaces_for_density

        script:
        """
        scil_transform_surface.py ${surfaces} ${affine_transfo}\
            ${sid}__surfaces_b0.vtk\
            --ants_warp ${warp_transfo}
        """
    }
}

surfaces_for_smooth
    .join(flow_mask_for_sf)
    .set{data_for_surface_smooth}

process D__Surface_Flow {
    cpus 1

    input:
    set sid, file(surf), file(mask)\
        from data_for_surface_smooth

    output:
    set sid, "${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.vtk"\
        into surface_for_intersections
    set sid, "${sid}__flow_${params.surf_flow_nb_step}_${params.surf_flow_step_size}.hdf5"\
        into surface_flow_lines

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
// *****************************************
// END COPY FROM SET-NF (with minor modification)
// *****************************************


surface_for_intersections
    .join(surface_flow_lines)
    .join(intersections_mask_for_sf)
    .join(surface_type_for_sf)
    .set{surface_and_map_for_intersections}

tractogram_for_intersections
    .combine(surface_and_map_for_intersections, by : 0)
    .set{data_for_intersections}


process F__Surface_Enhanced_Tractography {
    cpus 1

    input:
    set sid, file(tractogram), file(surf), file(flow), file(s_mask), file(s_type)\
        from data_for_intersections

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

// REMOVED SEED
process G__Concatenate_Intersection {
    cpus 1

    input:
    set sid, file(intersections)\
        from intersections_grouped_for_concatenate

    output:
    set sid, "${sid}__set_c_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.npz"\
        into intersections_for_connectivity, intersections_for_density

    script:
    command_lines=""
    intersections.each{
        command_lines +=" ${it} "
    }
    """
    scil_concatenate_surfaces_intersections.py ${command_lines}\
        --output_intersections ${sid}__set_c_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.npz
    """
}

surfaces_for_connectivity
    .join(labels_for_sf)
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
    set sid, "${sid}__set_connectivity_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.npy"\
        into matrices_for_stats

    script:
    """
    scil_surface_intersections_to_connectivity.py ${surf} ${intersections} ${labels}\
        ${sid}__set_connectivity_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.npy
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
    set sid, "${sid}__set_density_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.npy"\
        into density_for_stats

    script:
    """
    scil_surface_intersections_density.py ${surf} ${intersections}\
    ${sid}__set_density_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.npy --normalize_l1_to 1
    """
}

// streamlines_for_concatenate
//     .groupTuple()
//     .set{streamlines_grouped_for_concatenate}
//
// process Concatenate_Tractograms {
//     cpus 1
//
//     input:
//     set sid, file(tractograms) from streamlines_grouped_for_concatenate
//
//     output:
//     file "${sid}__set_${params.nb_seeds_per_surf}_c_filtered_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.fib"
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
//         ${sid}__set_${params.nb_seeds_per_surf}_c_filtered_${(int)params.minimum_length}_to_${(int)params.maximum_length}mm.fib
//     """
// }
