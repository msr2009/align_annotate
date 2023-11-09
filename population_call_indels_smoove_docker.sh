#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail #no POSIX

#requires docker, smoove docker image
#requires bcftools, bedtools, parallel

HELP(){
	echo "population_call_indels_smoove.sh"
	echo "wrapper script to call perform population calling with smoove" 
	echo "via docker. pipeline as suggested by brentp in smoove github."
	echo
	echo "syntax: population_call_indels_smoove.sh -g FASTA -b BAM -d DIR -t THREADS "
	echo "-g, --genome	path to FASTA genome"
	echo "--bamlist	 	list of bam files to be called. is input to parallel."
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
		--bamlist)
			BAM="$2"
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

NAME=`basename ${BAM%%.*}`

_FASTAFOLDER=`dirname ${GENOME}`
_GENOMEFASTA=`basename ${GENOME}`
SMOOVEDIR=${WORKINGDIR}/smoove/
SMOOVEVCF=${SMOOVEDIR}/${NAME}-smoove.genotyped.vcf.gz

#if smoove has been run before, there will be a .csi index 
#the presence of this causes smoove to stop running before 
#doing the bcftools splits below. So, let's remove everything 
#from the smoove directory before running

#smoove paste --name $cohort results-genotyped/*.vcf.gz

#do I have a useful gff file for worms



if [ -d ${SMOOVEDIR} ]; then
	echo "found existing smoove directory. overwriting."
	rm -f ${SMOOVEDIR}/*
fi

#want to parallelize this 

#smoove call --outdir results-smoove/ --exclude $bed --name $sample --fasta $reference_fasta -p 1 --genotype /path/to/$sample.bam
docker run \
               --platform linux/amd64 \
               --mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
               --mount type=bind,src=${WORKINGDIR},dst=/BAM \
               brentp/smoove smoove call \
               --name ${NAME} \
               --outdir /BAM/smoove/ \
               --fasta /FASTA/${_GENOMEFASTA} \
               --processes 1 \
			   --genotype \
               /BAM/${NAME}


#smoove merge --name merged -f $reference_fasta --outdir ./ results-smoove/*.genotyped.vcf.gz
docker run \
               --platform linux/amd64 \
               --mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
               --mount type=bind,src=${WORKINGDIR},dst=/BAM \
               brentp/smoove smoove merge \
               --name merged \
               --outdir /BAM/smoove/ \
               --fasta /FASTA/${_GENOMEFASTA} \
               --processes ${THREADS} \
               /BAM/smoove/*.genotyped.vcf.gz

#smoove genotype -d -x -p 1 --name $sample-joint --outdir results-genotped/ --fasta $reference_fasta --vcf merged.sites.vcf.gz /path/to/$sample.$bam
docker run \
               --platform linux/amd64 \
               --mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
               --mount type=bind,src=${WORKINGDIR},dst=/BAM \
               brentp/smoove smoove call \
               --name ${BAM} \
               --outdir /BAM/smoove/ \
               --fasta /FASTA/${_GENOMEFASTA} \
               --processes 1 \
			   --genotype \
               /BAM/${NAME}

#smoove annotate --gff Homo_sapiens.GRCh37.82.gff3.gz $cohort.smoove.square.vcf.gz | bgzip -c > $cohort.smoove.square.anno.vcf.gz
docker run \
               --platform linux/amd64 \
               --mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
               --mount type=bind,src=${WORKINGDIR},dst=/BAM \
               brentp/smoove smoove call \
               --name ${BAM} \
               --outdir /BAM/smoove/ \
               --fasta /FASTA/${_GENOMEFASTA} \
               --processes 1 \
			   --genotype \
               /BAM/${NAME}

