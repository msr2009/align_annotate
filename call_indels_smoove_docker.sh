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
SMOOVEVCF=${SMOOVEDIR}/${NAME}-smoove.genotyped.vcf.gz

#if smoove has been run before, there will be a .csi index 
#the presence of this causes smoove to stop running before 
#doing the bcftools splits below. So, let's remove everything 
#from the smoove directory before running

if [ -d ${SMOOVEDIR} ]; then
	echo "found existing smoove directory. overwriting."
	rm -f ${SMOOVEDIR}/*
fi

#run smoove within docker  
#(without duphold)
docker run \
               --platform linux/amd64 \
               --mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
               --mount type=bind,src=${WORKINGDIR},dst=/BAM \
               brentp/smoove smoove call --genotype\
               --name ${NAME} \
               --outdir /BAM/smoove/ \
               --fasta /FASTA/${_GENOMEFASTA} \
               --processes ${THREADS} \
               /BAM/${NAME}.srt.rmdup.bam

#run with duphold now
docker run \
               --platform linux/amd64 \
               --mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
               --mount type=bind,src=${WORKINGDIR},dst=/BAM \
               brentp/smoove smoove duphold \
               --vcf /BAM/smoove/${NAME}-smoove.genotyped.vcf.gz \
               --outvcf /BAM/smoove/${NAME}-smoove.genotyped.duphold.vcf.gz \
               --fasta /FASTA/${_GENOMEFASTA} \
               --processes ${THREADS} \
               /BAM/${NAME}.srt.rmdup.bam
