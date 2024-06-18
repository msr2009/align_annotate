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
ALIGNER="bwa"
BASECALLER="gatk"
BGVCF="none"
INDELS=0
TMPDIR="tmpdir"

#display help
HELP(){
	echo "align_annotate.sh"
	echo "an alignment pipeline (Matt Rich 2024)"
	echo
	echo "syntax: -d WORKING_DIRECTORY -1 READ1 -2 READ2 -g GENOME_FASTA -x PREFIX"
	echo
	echo "required arguments:"
	echo "-d, --dir	working directory"
	echo "-x, --prefix	prefix name for output"
	echo "-g, --genome	location of genome FASTA file"
	echo "-1		location of forward read FASTQ"
	echo "-2		location of reverse read FASTQ"
	echo
	echo "optional: alignment and basecalling:"
	echo "--aligner	alignment software to use 
		(bwa, bowtie2; default=bwa)"
	echo "--basecaller	basecalling software to use 
		(samtools, gatk; default=gatk)"
	echo "--call-indels	call indels with smoove and manta. 
		req's call_indels conda enviroment and python 2.7"
	echo
	echo "optional: filtration"	
	echo "--gq		genotype quality (GQ) threshold (default=15)"
	echo "--dp		read depth (DP) threshold (default =5)"
	echo "--bgvcf		multisample vcf with BGAF field for background filtering 
		must be made with 'make_background_vcf.py'"
	echo "--bgaf	background allele frequency threshold (default = 0.1)"
	echo "--strains_list	file containing all strains in population, one per line
		strains = folder names in directory."
	echo
	echo "optional: annotation"
	echo "-db	name of snpEff database for annotation (default=Caenorhabditis_elegans)"
	echo
	echo "optional: misc"
	echo "-t, --threads	number of threads to use for processes"
	echo "--tmpdir		directory for temporary files"
	echo "--params	parameters file. variables in file overwrite any defined here."
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
		--bgvcf)
			BGVCF="$2"
			shift 2
			;;
		--call-indels)
			INDELS=1
			shift 1
			;;
		--tmpdir)
			TMPDIR="$2"
			shift 2
			;;
		--gq)
			LOWQUAL="$2"
			shift 2
			;;
		--dp)
			LOWDEPTH="$2"
			shift 2
			;;
		--bgaf)
			BGAF="$2"
			shift 2
			;;
		--params_file)
			PARAMS="$2"
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

if [ -z ${TMPDIR} ]
then
		TMPDIR=${WORKING_DIR}/tempdir/
fi

#check architecture.
#intel procs can take advantage of hardware accel
#in gatk basecalling. 
ARCH=`uname -p`

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
echo TMPDIR: ${TMPDIR}
echo CPU: ${ARCH}
echo "######################################"

#make working directory if it doesn't exist
if [ ! -d ${WORKING_DIR} ]
then
		echo "Working directory (${WORKING_DIR}) does not exist. Creating."
		mkdir ${WORKING_DIR}	
fi

#make tempdir if it doesn't exist
if [ ! -d ${TMPDIR} ]
then
		echo "Temporary directory (${TMPDIR}) does not exist. Creating."
		mkdir ${TMPDIR}	
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
	if [ ${ARCH} = "x86_64" ]; then
		sh call_variants_gatk.sh -a ${_name}.srt.rmdup.bam -g ${GENOME} -t ${THREADS} 
	else
		sh call_variants_gatk_parallel.sh -a ${_name}.srt.rmdup.bam -g ${GENOME} -t ${THREADS} --parallel 500000
	fi
fi

if [ ${INDELS} = 1 ]
then
	echo
	echo "######################################"
	echo "CALLING INDELS"
	echo "######################################"

	sh indels_process_annotate.sh -d ${WORKING_DIR} -g ${GENOME} -x ${PREFIX} -t ${THREADS} --db ${DATABASE} --tmpdir ${TMPDIR}
fi
	
###filter vcfs
echo
echo "######################################"
echo "SOFT-FILTERING VCF"
echo "######################################"

if [ ${BGVCF} = "none" ]
then
	python soft-filter.py -v ${_name}.snp.vcf.gz
else
	python soft-filter.py -v ${_name}.snp.vcf.gz -b ${BGVCF}
fi

echo
echo "######################################"
echo "ANNOTATING SNPS VCF WITH SNPEFF"
echo "######################################"

#always call snp file
sh snpeff_annotation.sh --vcf ${_name}.snp.soft-filter.vcf.gz --db ${DATABASE}

