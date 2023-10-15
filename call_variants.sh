#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail

#REQUIRES: bcftools

HELP(){
        echo "call_variants.sh"
        echo "wrapper script for bcftools mpileup and call commands"
		echo 
		echo "required arguemnts:"
		echo "-a, --alignments	BAM or SAM file containing alignments"
		echo "-g, --genome	FASTA file for reference genome"
		echo "optional argument:"
		echo "--parallel	parallelize basecalling. argument = window size."
		echo "				REQUIRES: parallel, bedtools, faidx of genome."
		echo "-t, --threads	number of threads"
		echo "--gvcf		output GVCF instead of VCF"
		echo 
        echo "syntax: call_variants.sh [BAM/SAM] [FASTA]"
}

PARALLEL=0
TMPDIR="paralleltmpdir"
GENOME=""
THREADS=0
GVCF=0
ALIGNMENTS=""

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
				--parallel)
						PARALLEL="$2"
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
if [ ${PARALLEL} -eq 0 ]; then
	bcftools mpileup -Ou -f ${GENOME} ${ALIGNMENTS} | bcftools call -mv -f GQ -Oz -o ${_NAME}.snp.vcf.gz
else
	if [ -d ${TMPDIR} ]; then
			rm -rf ${TMPDIR}
	fi
	mkdir ${TMPDIR}
	bedtools makewindows -g ${GENOME}.fai -w ${PARALLEL} | bioawk -t '{print $1 ":" $2+1 "-" $3}' > ${TMPDIR}/tmp.windows
	parallel --jobs ${THREADS} bcftools mpileup -Ou -r {1} -f ${GENOME} ${ALIGNMENTS} '|' bcftools call -mv -f GQ -Oz -o ${TMPDIR}/tmp.{1}.vcf.gz :::: ${TMPDIR}/tmp.windows
	#then need to merge region'd vcfs
	bcftools concat --threads ${THREADS} -Oz ${TMPDIR}/tmp.*.vcf.gz | bcftools sort -Oz -o ${_NAME}.snp.vcf.gz
	
	rm -rf ${TMPDIR}
fi

#soft-filter heterozygous site (GT="hom") and 
#for low genotype quality ("GQ>=20")
#bcftools filter -s lowQual -i 'GQ>=10' ${_NAME}.snp.vcf.gz | bcftools filter -oZ -s Het -i 'GT=“hom"' -o ${_NAME}.snp.soft-filter.vcf.gz -Oz
