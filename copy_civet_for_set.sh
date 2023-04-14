#!/usr/bin/env bash
usage() { echo "$(basename $0) [-i civet] [-o output]" 1>&2; exit 1; }

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

CIVETINDIR=$(readlink -m ${i})
OUTDIR=$(readlink -m ${o})

echo "CIVET folder: ${CIVETINDIR}"
echo "Output folder: ${OUTDIR}"

echo "Building tree..."
cd ${CIVETINDIR}
for i in *;
do
    mkdir -p ${OUTDIR}/${i}

    # CIVET folders
    mkdir -p ${OUTDIR}/${i}/surfaces
    for cfile in gray_surface_left_81920.obj gray_surface_right_81920.obj white_surface_left_81920.obj white_surface_right_81920.obj;
    do
        FILE=`find ${CIVETINDIR}/${i} -name "*${cfile}" -print -quit`
        if [ -n "${FILE}" ]; then
            cp -L ${FILE} ${OUTDIR}/${i}/surfaces/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    mkdir -p ${OUTDIR}/${i}/transforms/linear
    FILE=`find ${CIVETINDIR}/${i} -name "*t1_tal.xfm" -print -quit`
    if [ -n "${FILE}" ]; then
        cp -L ${FILE} ${OUTDIR}/${i}/transforms/linear/
    else
        echo "WARNING! ${i} t1_tal.xfm was not found"
    fi

    mkdir -p ${OUTDIR}/${i}/segment
    FILE=`find ${CIVETINDIR}/${i} -name "*animal_labels.mnc" -print -quit`
    if [ -n "${FILE}" ]; then
        cp -L ${FILE} ${OUTDIR}/${i}/segment/
    else
        echo "WARNING! ${i} animal_labels.mnc was not found"
    fi

done
echo "Done"
