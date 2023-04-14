#!/usr/bin/env bash
usage() { echo "$(basename $0) [-i freesurfer] [-o output]" 1>&2; exit 1; }

while getopts "i:o:" args; do
    case "${args}" in
        i) i=${OPTARG};;
        o) o=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${i}" ] || [ -z "${o}" ]; then
    usage
fi

FSINDIR=$(readlink -m ${i})
OUTDIR=$(readlink -m ${o})

if [ -d "${FSINDIR}/results/" ]; then
    FSINDIR=${FSINDIR}/results
fi

echo "Freesurfer folder: ${FSINDIR}"
echo "Output folder: ${OUTDIR}"

echo "Building tree..."
cd ${FSINDIR}
for i in *;
do
    mkdir -p ${OUTDIR}/${i}

    # Freesurfer folders
    mkdir -p ${OUTDIR}/${i}/label
    for cfile in lh.aparc.annot rh.aparc.annot lh.aparc.a2009s.annot rh.aparc.a2009s.annot;
    do
        FILE=`find ${FSINDIR}/${i} -name ${cfile} -print -quit`
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/label/${cfile}
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/mri
    FILE=`find ${FSINDIR}/${i} -name "wmparc.mgz" -print -quit`
    if [ -n "${FILE}" ]; then
        cp -L ${FILE} ${OUTDIR}/${i}/mri/wmparc.mgz
    else
        FILE=`find ${FSINDIR}/${i} -name "wmparc.nii*" -print -quit`
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/mri/
        else
            echo "WARNING! ${i} wmparc was not found"
        fi
    fi

    mkdir -p ${OUTDIR}/${i}/surf
    for cfile in lh.pial lh.white rh.pial rh.white;
    do
        FILE=`find ${FSINDIR}/${i} -name ${cfile} -print -quit`
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/surf/${cfile}
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

done
echo "Done"
