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
		echo 
        echo "syntax: call_variants.sh [BAM/SAM] [FASTA]"
}

THREADS=1

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
	gatk HaplotypeCallerSpark -I ${ALIGNMENTS} -O ${_NAME}.vcf.gz -R ${GENOME} --spark-master local[${THREADS}]
else 
	gatk HaplotypeCaller -I ${ALIGNMENTS} -O ${_NAME}.vcf.gz -R ${GENOME}
fi

#after basecalling also create a vcf file containing only homozygous variants
#bcftools view -Oz -g hom ${_NAME}.vcf.gz > ${_NAME}.hom.vcf.gz

#after basecalling soft-filter variants for low quality (GQ>=20) and
#heterozygous genotypes
bcftools filter --threads ${THREADS} -s lowQual -i 'GQ>=20' ${_NAME}.vcf.gz | bcftools filter --threads ${THREADS} -s Het -i 'GT="hom"' -o ${_NAME}.soft-filter.vcf.gz -Oz
