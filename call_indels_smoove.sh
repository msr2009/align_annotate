#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail #no POSIX

#requires docker, smoove docker image
#requires bcftools, bedtools

HELP(){
	echo "call_indels_smoove.sh"
	echo "wrapper script to call smoove via docker"
	echo
	echo "syntax: call_indels_smoove.sh -g FASTA -n NAME -d DIR -t THREADS "
	echo "-g, --genome	path to FASTA genome"
	echo "-n, --name	name prefix for output (usually strain name)"
	echo "-d, --dir		path to working directory"
	echo "-t, --threads	number of threads to use for alignment"
	echo "-h, --help	show this help"
}

#default parameter values
THREADS=1

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

_FASTAFOLDER=`dirname ${GENOME}`
_GENOMEFASTA=`basename ${GENOME}`
SMOOVEDIR=${WORKINGDIR}/smoove/
SMOOVEVCF=${SMOOVEDIR}/${NAME}-smoove.genotyped.vcf

#if smoove has been run before, there will be a .csi index 
#the presence of this causes smoove to stop running before 
#doing the bcftools splits below. So, let's remove everything 
#from the smoove directory before running

if [ -d ${SMOOVEDIR} ]; then
	echo "found existing smoove directory. overwriting."
	rm -f ${SMOOVEDIR}/*
fi


##run smoove (without duphold)
#smoove call --genotype \
#			--name ${NAME} \
#			--outdir ${WORKINGDIR}/smoove/ \
#			--fasta ${_FASTAFOLDER}/${_GENOMEFASTA} \
#			--processes ${THREADS} \
#			${WORKINGDIR}/${NAME}.srt.rmdup.bam

#I've had trouble runnign the whole thing from smoove, but can do it in steps
#1) smoove call
smoove call -n ${NAME} -f ${GENOME} -p ${THREADS} -o ${SMOOVEDIR} ${WORKING_DIR}/${NAME}.srt.rmdup.bam
#2) svtyper-sso
gunzip ${SMOOVEVCF}.gz
svtyper-sso -i ${SMOOVEVCF} -o ${SMOOVEVCF%%.vcf}.genotyped.vcf --cores ${THREADS} -B ${WORKING_DIR}/${NAME}.srt.rmdup.bam
