#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail

#requires bowtie2, samtools

HELP(){
	echo "snpeff_annotation.sh"
	echo "wrapper script to annotate VCF file with snpEff"
	echo
	echo "syntax: snpeff_annotation.sh --db DATABASE --vcf VCF"
	echo "--db	snpEff database (e.g., 'Caenorhabditis_elegans')"
	echo "--vcf	VCF file"
	echo "-h, --help	show this help"
}

#PARSE OPTIONS
while [ $# -gt 0 ]; do
	case "$1" in
		--db)
			DATABASE="$2"
			shift 2
			;;
		--vcf)
			VCF="$2"
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

ANN_NAME=${VCF}.snpeff
if [[ ${VCF} == *.vcf.gz ]]; then
	ANN_NAME=${VCF%%.vcf.gz}.snpeff
elif [[ ${VCF} == *.vcf ]]; then
	ANN_NAME=${VCF%%.vcf}.snpeff
fi



###perform snpeff annotation
snpEff eff -noStats -ud 0 ${DATABASE} ${VCF} > ${ANN_NAME}.vcf
###convert annotated VCF to a more human-readable format
python snpeff2tsv.py --vcf ${ANN_NAME}.vcf > ${ANN_NAME}.tsv

