#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, samtools, blast, bioawk

# Default values for blank parameters
GENOME=
PREFIX=
THREADS=1
DATABASE="Caenorhabditis_elegans"
ALIGNER="bwa"
BASECALLER="gatk"
SMOOVE=1
SNPEFF_INPUT=""

#display help
HELP(){
	echo "align_annotate.sh"
	echo "an alignment pipeline (Matt Rich 2021)"
	echo
	echo "syntax: -d WORKING_DIRECTORY -1 READ1 -2 READ2 -g GENOME_FASTA -x PREFIX"
	echo "required options:"
	echo "-d, --dir	working directory"
	echo "-x, --prefix	prefix name for output"
	echo "-g, --genome	location of genome FASTA file"
	echo "-1		location of forward read FASTQ"
	echo "-2		location of reverse read FASTQ"
	echo
	echo "optional parameters:"
	echo "--aligner		alignment software to use" 
	echo "					(bwa, bowtie2; default=bwa)"
	echo "--basecaller	basecalling software to use" 
	echo "					(samtools, gatk; default=gatk)"
	echo "-db	name of snpEff database for annotation (default=Caenorhabditis_elegans)"
	echo "-t, --threads	number of threads to use for processes"
	echo "--no-smoove	do not use smoove to call indels"
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
		-1)
            READ1="$2"
            shift 2
            ;;
        -2)
            READ2="$2"
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
		--aligner)
			ALIGNER="$2"
			shift 2
			;;
		--basecaller)
			BASECALLER="$2"
			shift 2
			;;
		--no-smoove)
			SMOOVE=0
			shift 1
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

if [ -z $READ1 ]
then
	echo "ERROR: must set read1 path with -1"
	echo
	HELP
	exit
fi

if [ -z $READ2 ]
then
	echo "ERROR: must set read2 path with -2"
	echo
	HELP
	exit
fi

##############
###MAIN SCRIPT	
##############
_name=${WORKING_DIR}/${PREFIX}

echo
echo "######################################"
echo NAME: ${PREFIX}
echo GENOME: ${GENOME}
echo READ1: ${READ1}
echo READ2: ${READ2}
echo WORKING DIRECTORY: ${WORKING_DIR}
echo ALIGNER: ${ALIGNER}
echo BASECALLER: ${BASECALLER}
echo SNPEFF DATABASE: ${DATABASE}
echo "######################################"

#make working directory if it doesn't exist
if [ ! -d ${WORKING_DIR} ]
then
		echo "Working directory (${WORKING_DIR}) does not exist. Creating."
		mkdir ${WORKING_DIR}	
fi

echo
echo "######################################"
echo "ALIGNING READS TO GENOME"
echo "######################################"

if [ ${ALIGNER} = "bwa" ]
then
	sh aln_bwa.sh -g ${GENOME} -x ${_name} -t ${THREADS} -1 ${READ1} -2 ${READ2}
elif [ ${ALIGNER} = "bowtie2" ]
then
	sh aln_bt2.sh -g ${GENOME} -x ${_name} -t ${THREADS} -1 ${READ1} -2 ${READ2}
fi

echo
echo "######################################"
echo "PROCESSING ALIGNED READS"
echo "######################################"

sh process_alignment.sh -i ${_name}.bam -t ${THREADS}

echo
echo "######################################"
echo "BASE CALLING"
echo "######################################"

if [ ${BASECALLER} = "samtools" ]
then
	echo "basecalling with samtools"
	sh call_variants.sh ${_name}.srt.rmdup.bam ${GENOME} -t ${THREADS}
elif [ ${BASECALLER} = "gatk" ]
then
	echo "basecalling with gatk"
	sh call_variants_gatk.sh -a ${_name}.srt.rmdup.bam -g ${GENOME} -t ${THREADS}
fi


###HARD CODING THIS IN FOR NOW BECAUSE MY NO-SMOOVE OPTION ISN'T WORKING....
#SMOOVE=0

if [ ${SMOOVE} = 1 ]
then
	echo
	echo "######################################"
	echo "CALLING INDELS"
	echo "######################################"

	sh call_indels_smoove.sh -d ${WORKING_DIR}/smoove/ -n ${PREFIX} -g ${GENOME} -t ${THREADS}
	
	#concatenate snp and indel calls into same vcf file
	bcftools concat -a -Oz -o ${_name}.all.soft-filter.vcf.gz ${_name}.snv.soft-filter.vcf.gz ${_name}.dup.vcf.gz ${_name}.del.vcf.gz
	SNPEFF_INPUT=${_name}.all.soft-filter.vcf.gz
else
	echo
	echo "######################################"
	echo "SKIPPING INDEL CALLING STEP"
	echo "######################################"
	SNPEFF_INPUT=${_name}.snp.soft-filter.vcf.gz
fi

echo
echo "######################################"
echo "ANNOTATING VCF WITH SNPEFF"
echo "######################################"

sh snpeff_annotation.sh --vcf ${SNPEFF_INPUT} --db ${DATABASE}


