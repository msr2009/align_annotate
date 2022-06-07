#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

HELP(){
        echo "call_variants.sh"
        echo "wrapper script for GATK HaplotypeCaller"
        echo "required arguemnts:"
		echo "-a, --alignments	BAM or SAM file containing alignments"
		echo "-g, --genome	FASTA file for reference genome"
		echo "optional argument:"
		echo "-t, --threads	number of threads"
		echo "--gvcf		output GVCF instead of VCF"
		echo 
        echo "syntax: call_variants_gatk.sh [BAM/SAM] [FASTA]"
}

THREADS=1
GVCF=0

while [ $# -gt 0 ]; do
        case "$1" in
                -a|--alignments)
                        ALIGNMENTS="$2"
                        shift 2
                        ;;
                -g|--genome)
                        GENOME="$2"
                        shift 2
                        ;;
                -t|--threads)
                        THREADS="$2"
                        shift 2
                        ;;
				--gvcf)
						GVCF=1
						shift 1
						;;
                -h|--help) #HELP!
                        HELP
                        exit;;
                *)
                        HELP
                        exit;;
        esac
done

_NAME="${ALIGNMENTS%%.*}"

#then perform basecalling 
echo "alignments: ${ALIGNMENTS}"
echo "genome: ${GENOME}"

if [ $THREADS -ne 1 ]; then
	if [ $GVCF -eq 0 ]; then	
		gatk HaplotypeCallerSpark -I ${ALIGNMENTS} -O ${_NAME}.snp.vcf.gz -R ${GENOME} --spark-master local[${THREADS}]
	else
		gatk HaplotypeCallerSpark -I ${ALIGNMENTS} -O ${_NAME}.g.vcf.gz -R ${GENOME} --spark-master local[${THREADS}] -ERC GVCF
	fi
else 
	if [ $GVCF -eq 0 ]; then
		gatk HaplotypeCaller -I ${ALIGNMENTS} -O ${_NAME}.snp.vcf.gz -R ${GENOME}
	else
		gatk HaplotypeCaller -I ${ALIGNMENTS} -O ${_NAME}.g.vcf.gz -R ${GENOME} -ERC GVCF
	fi
fi

#after basecalling soft-filter variants for low quality (GQ>=20) and
#heterozygous genotypes
#only do this if not outputting gVCF
if [ $GVCF -eq 0 ]; then
	bcftools filter --threads ${THREADS} -s lowQual -i 'GQ>=20' ${_NAME}.snp.vcf.gz | bcftools filter --threads ${THREADS} -s Het -i 'GT="hom"' -o ${_NAME}.snp.soft-filter.vcf.gz -Oz
	bcftools index ${_NAME}.snp.soft-filter.vcf.gz
fi
