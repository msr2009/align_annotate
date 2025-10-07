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
	echo "for performing indel analysis and annotation of aligned reads" 
	echo "calls indels with manta and smoove"
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

###there's some weirdness here:
###smoove and manta run via python2 (this is one reason by smoove is usually in a docker image)
###pysam needs python3
###so, we need to swap conda environments in the middle of this script. 

###I CAN USE `conda run` for this! much nicer!
###I also have to change "call_indels_smoove.sh" to not use docker anymore
###(make "call_indels_smoove_docker" and "call_indels_smoove")

#conda run -n call_indels  call_indels_smoove.sh -d ${WORKING_DIR} -n ${PREFIX} -g ${GENOME} -t ${THREADS}
#conda run -n call_indels  call_indels_manta.sh -d ${WORKING_DIR} -n ${PREFIX} -g ${GENOME} -t ${THREADS}


#we're going to assume we're starting in align_annotate, and swap first into call_indels
echo "swapping into python2 environment (call_indels) for indel calling"
eval "$(conda shell.bash hook)"
conda activate call_indels

#call indels with smoove
sh call_indels_smoove.sh -d ${WORKING_DIR} -n ${PREFIX} -g ${GENOME} -t ${THREADS}
#call indels with manta
sh call_indels_manta.sh -d ${WORKING_DIR} -n ${PREFIX} -g ${GENOME} -t ${THREADS}

###here's where we swap back to align_annotate
echo "moving back to align_annotate environment"
#eval "$(conda shell.bash hook)"
conda activate align_annotate

#concatenate smoove and manta indels
bcftools concat -a -o ${_name}.allSV.vcf.gz -Oz ${WORKING_DIR}/smoove/${PREFIX}-smoove.genotyped.duphold.vcf.gz ${WORKING_DIR}/manta/results/variants/diploidSV.vcf.gz
bcftools index ${_name}.allSV.vcf.gz

#filter using soft-filter.py
python soft-filter.py -v ${_name}.allSV.vcf.gz

#extract files for dels, dups, and insertions
bcftools filter -i 'INFO/SVTYPE="DEL" & INFO/SVLEN>-10000' -o ${_name}.del.soft-filter.vcf.gz -Oz ${_name}.allSV.soft-filter.vcf.gz
bcftools filter -i 'INFO/SVTYPE="DUP" & INFO/SVLEN<10000' -o ${_name}.dup.soft-filter.vcf.gz -Oz ${_name}.allSV.soft-filter.vcf.gz
bcftools filter -i 'INFO/SVTYPE="INS" & INFO/SVLEN<10000' -o ${_name}.ins.soft-filter.vcf.gz -Oz ${_name}.allSV.soft-filter.vcf.gz

echo
echo "######################################"
echo "ANNOTATING SV VCFS WITH SNPEFF"
echo "######################################"

#annotate dup, del, and ins files
sh snpeff_annotation.sh --vcf ${_name}.dup.soft-filter.vcf.gz --db ${DATABASE}
sh snpeff_annotation.sh --vcf ${_name}.del.soft-filter.vcf.gz --db ${DATABASE}
sh snpeff_annotation.sh --vcf ${_name}.ins.soft-filter.vcf.gz --db ${DATABASE}
