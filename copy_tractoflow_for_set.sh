#!/usr/bin/env bash
usage() { echo "$(basename $0) [-i tractoflow] [-o output]" 1>&2; exit 1; }

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

TFINDIR=$(readlink -m ${i})
OUTDIR=$(readlink -m ${o})

if [ -d "${TFINDIR}/results/" ]; then
    TFINDIR=${TFINDIR}/results
fi

echo "Tractoflow/results folder: ${TFINDIR}"
echo "Output folder: ${OUTDIR}"

echo "Building tree..."
cd ${TFINDIR}
for i in *;
do
    mkdir -p ${OUTDIR}/${i}

    # Tractoflow folders
    mkdir -p ${OUTDIR}/${i}/DTI_Metrics
    for cfile in fa.nii.gz;
    do
        FILE=`find ${TFINDIR}/${i} -name "*__${cfile}" -print `
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/DTI_Metrics/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/FODF_Metrics
    for cfile in fodf.nii.gz;
    do
        FILE=`find ${TFINDIR}/${i} -name "*__${cfile}" -print `
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/FODF_Metrics/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/PFT_Maps
    for cfile in map_exclude.nii.gz map_include.nii.gz;
    do
        FILE=`find ${TFINDIR}/${i} -name "*__${cfile}" -print `
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/PFT_Maps/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/Register_T1
    for cfile in output0GenericAffine.mat output1InverseWarp.nii.gz;
    do
        FILE=`find ${TFINDIR}/${i} -name "*__${cfile}" -print `
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/Register_T1/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

done
echo "Done"
