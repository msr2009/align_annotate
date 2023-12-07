#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail #no POSIX

#requires manta, bedtools

HELP(){
	echo "call_indels_manta.sh"
	echo "wrapper script to call manta"
	echo
	echo "syntax: call_indels_manta.sh -g FASTA -n NAME -d DIR -t THREADS"
	echo "-g, --genome	path to FASTA genome"
	echo "-n, --name	name prefix for output (usually strain name)"
	echo "-d, --dir	path to working directory"
	echo "-p, --path-to-manta	path to manta bin (auto-detected as `which configManta.py`"
	echo "-t, --threads	number of threads to use for alignment"
	echo "-h, --help	show this help"
}

#default parameter values
THREADS=1
MANTAPATH=`which configManta.py`

#PARSE OPTIONS
while [ $# -gt 0 ]; do
	case "$1" in
		-g|--genome)
			GENOME="$2"
			shift 2
			;;
		-n|--name)
			NAME="$2"
			shift 2
			;;
		--path-to-manta)
			MANTAPATH="$2"
			shift 2
			;;
		-d|--dir) 
			WORKINGDIR="$2"
			shift 2
			;;
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-h|--help) #HELP!
			HELP
			exit;;
		*)
			HELP
			exit;;
	esac
done

MANTADIR=${WORKINGDIR}/manta/

if [ -d ${MANTADIR} ]; then
	echo "found existing manta directory. overwriting."
	rm -f ${MANTADIR}/*
fi

python ${MANTAPATH} \
	--bam ${WORKINGDIR}/${NAME}.srt.rmdup.bam \
	--referenceFasta ${GENOME} \
	--runDir ${MANTADIR}

python ${MANTADIR}/runWorkflow.py -j ${THREADS}

###should probably add some post-processing to manta vcf...change to <DUP/DEL> etc...
