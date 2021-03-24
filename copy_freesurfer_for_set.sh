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
        if [ -f ${FSINDIR}/${i}/label/${cfile} ]; then
            cp -L ${FSINDIR}/${i}/label/${cfile}  ${OUTDIR}/${i}/label/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/mri
    for cfile in wmparc;
    do
        if [ -f ${FSINDIR}/${i}/mri/${cfile}* ]; then
            cp -L ${FSINDIR}/${i}/mri/${cfile}*  ${OUTDIR}/${i}/mri/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/surf
    for cfile in lh.pial lh.white rh.pial rh.white;
    do
        if [ -f ${FSINDIR}/${i}/surf/${cfile} ]; then
            cp -L ${FSINDIR}/${i}/surf/${cfile}  ${OUTDIR}/${i}/surf/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

done
echo "Done"
