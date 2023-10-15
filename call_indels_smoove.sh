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


#run smoove within docker 
docker run \
		--platform linux/amd64 \
		--mount type=bind,src=${_FASTAFOLDER},dst=/FASTA \
	  	--mount type=bind,src=${WORKINGDIR},dst=/BAM \
		brentp/smoove smoove call -d --genotype \
		--name ${NAME} \
		--outdir /BAM/smoove/ \
		--fasta /FASTA/${_GENOMEFASTA} \
		/BAM/${NAME}.srt.rmdup.bam


### THIS iS CURRENTLY ALL WRONG
### MAKE THREE FILES:
#1) DELs soft-filtered for DHBFC<0.25
#2) DUPs soft-filtered for DHBFC>1.75
#3) everything else (not DEL or DUP)
#then we can concat everything together with the snvs vcf
SMOOVEVCF=${WORKINGDIR}/smoove/${NAME}-smoove.genotyped.vcf.gz

bcftools filter -i 'SVTYPE="DEL" & DHBFC<0.25 & SVLEN>-50000 & SVLEN<50000' -Oz -o ${WORKINGDIR}/${NAME}.del.vcf.gz ${SMOOVEVCF}
bcftools index ${WORKINGDIR}/${NAME}.del.vcf.gz

bcftools filter -i 'SVTYPE="DUP" & DHBFC>1.75 & SVLEN>-50000 & SVLEN<50000' -Oz -o ${WORKINGDIR}/${NAME}.dup.vcf.gz ${SMOOVEVCF}
bcftools index ${WORKINGDIR}/${NAME}.dup.vcf.gz

bcftools filter -i 'SVLEN>=50000 | SVLEN<=-50000' -Oz -o ${WORKINGDIR}/${NAME}.longSV.vcf.gz ${SMOOVEVCF}
