#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail #no POSIX

#requires smoove, parallel

HELP(){
	echo "population_call_indels_smoove.sh"
	echo "wrapper script to call perform population calling with smoove" 
	echo "via docker. pipeline as suggested by brentp in smoove github."
	echo
	echo "syntax: population_call_indels_smoove.sh -g FASTA -b BAM -d DIR -t THREADS "
	echo "-g, --genome	path to FASTA genome"
	echo "--bam-list	 	list of bam files to be called. is input to parallel."
	echo "-d, --dir		path to working directory"
	echo "-t, --threads	number of threads to use for alignment"
	echo "-o, --out-name	name output vcf"
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
		--bam-list)
			BAMLIST="$2"
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
		-o|--output-name)
			OUTNAME="$2"
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

OUTVCF=${OUTNAME%%.*}

#smoove call --outdir results-smoove/ --exclude $bed --name $sample --fasta $reference_fasta -p 1 --genotype /path/to/$sample.bam
#use parallel to parallelize across all threads. call runs best with only one
parallel -a ${BAMLIST} echo smoove call --outdir ${WORKINGDIR}/smoove --name {/.} --fasta ${GENOME} -p 1 --genotype {}

#smoove merge --name merged -f $reference_fasta --outdir ./ results-smoove/*.genotyped.vcf.gz
smoove merge --name merged -f ${GENOME} --outdir ${WORKINGDIR}/smoove -p ${THREADS} ${WORKINGDIR}/smoove/*.genotyped.vcf.gz

#smoove genotype -d -x -p 1 --name $sample-joint --outdir results-genotped/ --fasta $reference_fasta --vcf merged.sites.vcf.gz /path/to/$sample.$bam
for bam in `cat ${BAMLIST}`; do
	echo "genotyping ${bam} at sites in ${WORKINGDIR}/smoove/merged.sites.vcf.gz"
	NAME=`basename ${bam%%.*}`
	smoove genotype -d -x -p ${THREADS} --name ${NAME}-joint --outdir ${WORKINGDIR}/smoove/genotyped/ --fasta ${GENOME} --vcf ${WORKINGDIR}/smoove/merged.sites.vcf.gz ${bam}
done

#smoove paste --name $cohort results-genotyped/*.vcf.gz
smoove paste --name ${OUTVCF} ${WORKINGDIR}/smoove/genotyped/*.vcf.gz

#smoove annotate --gff Homo_sapiens.GRCh37.82.gff3.gz $cohort.smoove.square.vcf.gz | bgzip -c > $cohort.smoove.square.anno.vcf.gz
