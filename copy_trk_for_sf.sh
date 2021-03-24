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

echo "Tractoflow/results folder: ${TFINDIR}"
echo "Output folder: ${OUTDIR}"

echo "Building tree..."
cd ${TFINDIR}
for i in *;
do
    mkdir -p ${OUTDIR}/${i}

    # Tractoflow folders
    cp -L ${TFINDIR}/${i}/*/*.trk  ${OUTDIR}/${i}/


done
echo "Done"
