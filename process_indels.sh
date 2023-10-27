#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail #no POSIX

#requires docker, smoove docker image
#requires bcftools, bedtools

HELP(){
	echo "prcoess_indels.sh"
	echo "this script filters smoove output to generate vcf files"
	echo "containing high-quality SV calls"
	echo
	echo "syntax: process_indels.sh  "
	echo "--vcf		smoove vcf (usually NAME-smoove.genotyped.vcf.gz)"
	echo "--d, --dir	working directory"
	echo "--sv_length	maximum length of SV to add to dup/del files."
	echo "		longer SVs will be put in longSV file. default=50000"
	echo "-t, --threads	number of threads to use for alignment"
	echo "-h, --help	show this help"
}

#default parameter values
THREADS=1
SV_LENGTH=50000


#PARSE OPTIONS
while [ $# -gt 0 ]; do
	case "$1" in
		--vcf)
			SMOOVEVCF="$2"
			shift 2
			;;
		--sv_length)
			SV_LENGTH="$2"
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

NAME=`basename ${SMOOVEVCF} | cut -d "-" -f 1`

bcftools filter -i "SVTYPE='DEL' & DHFFC<0.3 & SVLEN>-${SV_LENGTH} & SVLEN<${SV_LENGTH}" -Oz -o ${WORKINGDIR}/${NAME}.del.vcf.gz ${SMOOVEVCF}
bcftools index ${WORKINGDIR}/${NAME}.del.vcf.gz

bcftools filter -i "SVTYPE='DUP' & DHBFC>1.7 & SVLEN>-${SV_LENGTH} & SVLEN<${SV_LENGTH}" -Oz -o ${WORKINGDIR}/${NAME}.dup.vcf.gz ${SMOOVEVCF}
bcftools index ${WORKINGDIR}/${NAME}.dup.vcf.gz

bcftools filter -i "SVLEN>=${SV_LENGTH} | SVLEN<=-${SV_LENGTH}" -Oz -o ${WORKINGDIR}/${NAME}.longSV.vcf.gz ${SMOOVEVCF}
bcftools index ${WORKINGDIR}/${NAME}.longSV.vcf.gz ${SMOOVEVCF}
