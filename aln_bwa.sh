#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, samtools

HELP(){
	echo "aln_bwa.sh"
	echo "wrapper script to align array reads to a genome using bwa-mem"
	echo "does _not_ perform post-alignment processsing"
	echo
	echo "syntax: aln_bwa.sh -g GENOME -x PREFIX -1 READ1 -2 READ2 -t THREADS"
	echo "-g, --genome	location of genome FASTA"
	echo "-x, --prefix	prefix for naming output (should contain full path)"
	echo "-1	location of forward read FASTQ"
	echo "-2	location of reverse read FASTQ"
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
		-x|--prefix)
			PREFIX="$2"
			shift 2
			;;
		-1) 
			READ1="$2"
			shift 2
			;;
		-2)
			READ2="$2"
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

RGID=`basename ${PREFIX}`

#check if genome has been indexed; in not, do it.
if [ ! -f ${GENOME}.bwt ]; then
		echo "creating bwa index of ${GENOME}"
		bwa index ${GENOME}
fi

bwa mem -t ${THREADS} -R "@RG\tID:${RGID}\tLP:lib\tPL:ILLUMINA\tPU:barcode\tSM:${RGID}" ${GENOME} ${READ1} ${READ2} | samtools view -bS - > ${PREFIX}.bam

