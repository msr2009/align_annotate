#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, samtools

HELP(){
	echo "aln_bt2.sh"
	echo "wrapper script to align array reads to a genome"
	echo "does _not_ perform post-alignment processsing"
	echo
	echo "syntax: aln_bt2.sh -g GENOME -x PREFIX -1 READ1 -2 READ2"
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

#create bowtie2 index
RGID=`basename ${PREFIX}`


IDX= ${GENOME%.*}
if [ ! -f ${IDX}.1.bt2 ]
		echo "building bowtie2 index of ${GENOME}"
		bowtie2-build --threads ${THREADS} ${GENOME} ${IDX}
fi

###align reads and convert to bam
bowtie2 -t -p ${THREADS} -x ${IDX} -1 ${READ1} -2 ${READ2} --rg-id ${RGID} --rg LP:lb --rg PL:ILLUMINA --rg PU:barcode --rg SM:${RGID} | samtools view -bS > ${PREFIX}.bam

