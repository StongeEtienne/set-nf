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
    for cfile in gray_surface_left_81920.obj gray_surface_right_81920.obj white_surface_left_81920.obj white_surface_right_81920.obj;
    do
        if [ -f ${CIVETINDIR}/${i}/surfaces/*${cfile} ]; then
            cp -L ${CIVETINDIR}/${i}/surfaces/*${cfile}  ${OUTDIR}/${i}/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    for cfile in t1_tal.xfm;
    do
        if [ -f ${CIVETINDIR}/${i}/transforms/linear/*${cfile} ]; then
            cp -L ${CIVETINDIR}/${i}/transforms/linear/*${cfile}  ${OUTDIR}/${i}/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

    for cfile in animal_labels.mnc;
    do
        if [ -f ${CIVETINDIR}/${i}/segment/*${cfile} ]; then
            cp -L ${CIVETINDIR}/${i}/segment/*${cfile}  ${OUTDIR}/${i}/
        else
            echo "WARNING! ${i} ${cfile} was not found"
        fi
    done

done
echo "Done"
