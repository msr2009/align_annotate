#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, bwa, samtools
#requires snpeff, snpeff2tsv.py

# Default values for blank parameters
GENOME=
FASTQ_DIR=
OUTPUT_DIR=
INFOFILE=
OVERWRITE=0
THREADS=1
DATABASE="Caenorhabditis_elegans"
ALIGNER="bwa"
BASECALLER="gatk"
WORKING_DIR=

#display help
HELP(){
	echo "batch_genome_alignment.sh"
	echo "a wrapper to an  alignment pipeline (Matt Rich 2021)"
	echo
	echo "takes tab-delimited file containing strain information"	
	echo "and fastq names and calls genome_alignment.sh for each"
	echo "strain."
	echo
	echo "syntax: --fastq_dir FASTQDIR --output_dir OUTPUTDIR --map INFOFILE --genome GENOME"
	echo "required options:"
	echo "-f, --fastq_dir 	directory containing fastq files from sequencing run"
	echo "-o, --output_dir	directory where each strain output file be added"
	echo "-g, --genome		location of genome FASTA file"
	echo "-m, --map		file containing strain and fastq information"	
	echo
	echo "optional parameters:"
	echo "--aligner		alignment software to use" 
	echo "				(bwa, bowtie2; default=bwa)"
	echo "--basecaller		basecalling software to use" 
	echo "				(samtools, gatk; default=gatk)"
	echo "-db			name of snpEff database for annotation" 
	echo "				(default=Caenorhabditis_elegans)"
	echo "-t, --threads		number of threads to use for processes (default=1)"
	echo "--overwrite		overwrite existing output if it already exists"
	echo "-h, --help		show this help"
}

if [ $# -eq 0 ]
then
	HELP
	exit
fi	

#PARSE OPTIONS
while [ $# -gt 0 ]; do
    case "$1" in
		-f|--fastq_dir)
			FASTQ_DIR="$2"
			shift 2
			;;
		-o|--output_dir)
	    	OUTPUT_DIR="$2"
	    	shift 2
	    	;;
		-m|--map)
			INFOFILE="$2"
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
		--overwrite)
			OVERWRITE=1
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
if [ -z $FASTQ_DIR ]
then
	echo "ERROR: must set fastq directory path with -f"
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

if [ -z "$INFOFILE" ]
then
	echo "ERROR: must set info file with -m"
	echo 
	HELP
	exit
fi	

if [ -z $OUTPUT_DIR ]
then
	echo "ERROR: must set output directory path with -o"
	echo
	HELP
	exit
fi	

##############

##check if output directory exists; if not make it.
if [ ! -d ${OUTPUT_DIR} ]
then
		echo "Output directory (${OUTPUT_DIR}) does not exist. Creating it."
		mkdir ${OUTPUT_DIR}
fi

#NOTE: this step relies on static positions for columns in
#the mapping file from the HCI sequence core. 

for strain in `awk '{FS="\t";OFS=","}NR>1{print $18, $19}' "${INFOFILE}"`; do
		
		FASTQ_PREFIX=$(echo $strain | cut -d "," -f 1)
		SAMPLE_NAME=$(echo $strain | cut -d "," -f 2 | sed 's/#//g')
		WORKING_DIR=${OUTPUT_DIR}/${SAMPLE_NAME}/

		if [ ${OVERWRITE} -eq 0 ] && [ -f ${WORKING_DIR}/${SAMPLE_NAME}.snpeff.tsv ]
		then
			echo "########################################################"
			echo "It looks like ${SAMPLE_NAME} has already been processed."
			echo "I'm going to skip it"	
			echo "########################################################"
			continue
		fi

		READ1=$(find ${FASTQ_DIR}/${FASTQ_PREFIX}_*R1_001.fastq.gz)
		READ2=$(find ${FASTQ_DIR}/${FASTQ_PREFIX}_*R2_001.fastq.gz)

		sh align_annotate.sh -d ${WORKING_DIR} -x ${SAMPLE_NAME} -g ${GENOME} \
								-1 ${READ1} -2 ${READ2} -t ${THREADS} \
								--aligner ${ALIGNER} --basecaller ${BASECALLER}

		#after processing, move FASTQ files into output folder
		mv ${READ1} ${WORKING_DIR}/
		mv ${READ2} ${WORKING_DIR}/
			
done
