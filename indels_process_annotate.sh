#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail #no POSIX

#requires bowtie2, samtools, blast, bioawk

# Default values for blank parameters
GENOME=
PREFIX=
THREADS=1
DATABASE="Caenorhabditis_elegans"
SNPEFF_INPUT=""
TMPDIR=

#display help
HELP(){
	echo "indel_process_annotate.sh"
	echo "for performing smoove analysis and annotation of aligned reads" 
	echo "(Matt Rich 2023)"
	echo
	echo "syntax: -d WORKING_DIRECTORY -1 READ1 -2 READ2 -g GENOME_FASTA -x PREFIX"
	echo "required options:"
	echo "-d, --dir	working directory"
	echo "-x, --prefix	prefix name for output"
	echo "-g, --genome	location of genome FASTA file"
	echo
	echo "optional parameters:"
	echo "-db	name of snpEff database for annotation (default=Caenorhabditis_elegans)"
	echo "-t, --threads	number of threads to use for processes"
	echo "--tmpdir	directory for temporary files"
	echo "-h, --help	show this help"
}

if [ $# -eq 0 ]
then
	HELP
	exit
fi	

#PARSE OPTIONS
while [ $# -gt 0 ]; do
    case "$1" in
		-d|--dir)
			WORKING_DIR="$2"
			shift 2
			;;
		-x|--prefix)
	    	PREFIX="$2"
	    	shift 2
	    	;;
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-g|--genome)
			GENOME="$2"
			shift 2
			;;
		--db)
			DATABASE="$2"
			shift 2
			;;
		--tmpdir)
			TMPDIR="$2"
			shift 2
			;;
		-h|--help) #HELP ME PLEASE!
			HELP
			exit;;
		*)
			HELP
            exit;;
    esac
done

#check that required parameters are set
if [ -z $WORKING_DIR ]
then
	echo "ERROR: must set working directory path with -d"
	echo
	HELP
	exit
fi	

if [ -z $GENOME ]
then
	echo "ERROR: must set path to genome with -g"
	echo
	HELP
	exit
fi	

if [ -z $PREFIX ]
then
	echo "ERROR: must set prefix with -x"
	echo 
	HELP
	exit
fi	

##############
###MAIN SCRIPT	
##############
_name=${WORKING_DIR}/${PREFIX}

if [ -z ${TMPDIR} ]
then
		TMPDIR=${WORKING_DIR}/tempdir/
fi

if [ ! -d ${TMPDIR} ]
then
		echo "Temporary directory (${TMPDIR}) does not exist. Creating."
		mkdir ${TMPDIR}	
fi

echo
echo "######################################"	
echo "CALLING INDELS"
echo "######################################"

#call indels with smoove
sh call_indels_smoove.sh -d ${WORKING_DIR} -n ${PREFIX} -g ${GENOME} -t ${THREADS}
	
#process smoove vcf into del and dup files
sh process_indels.sh -d ${WORKING_DIR} --vcf ${WORKING_DIR}/smoove/${PREFIX}-smoove.genotyped.duphold.vcf.gz

echo
echo "######################################"
echo "ANNOTATING VCF WITH SNPEFF"
echo "######################################"

#always call snp file
sh snpeff_annotation.sh --vcf ${_name}.snp.soft-filter.vcf.gz --db ${DATABASE}

#if dup and del files exist, also call those
if [ -f ${_name}.dup.vcf.gz ]; then
	sh snpeff_annotation.sh --vcf ${_name}.dup.vcf.gz --db ${DATABASE}
fi

if [ -f ${_name}.del.vcf.gz ]; then
	sh snpeff_annotation.sh --vcf ${_name}.del.vcf.gz --db ${DATABASE}
fi
