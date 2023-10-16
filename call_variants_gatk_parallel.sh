#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail

HELP(){
        echo "call_variants.sh"
        echo "wrapper script for GATK HaplotypeCaller"
        echo "required arguemnts:"
		echo "-a, --alignments	BAM or SAM file containing alignments"
		echo "-g, --genome	FASTA file for reference genome"
		echo "optional argument:"
		echo "--parallel	parallelize basecalling with gnu parallel"
		echo "		argument: window size for chunking genome"
		echo "			REQUIRES: parallel, bedtools, genome faidx"
		echo "-t, --threads	number of threads"
		echo "--gvcf		output GVCF instead of VCF"
		echo "--tmpdir	directory for temporary files"
		echo 
        echo "syntax: call_variants_gatk.sh [BAM/SAM] [FASTA]"
}

THREADS=1
GVCF=0
_TMPDIR=
PARALLEL=0

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
				--tmpdir)
						_TMPDIR="$2"
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

if [ -z $_TMPDIR ]
then
		_TMPDIR=`dirname ${_NAME}`/tempdir/
fi

##for some reason parallel stopped working 
##getting a GATK access denied error when it tries
##to write out... so set PARALLEL to 0. zannen.
#PARALLEL=0
#setting threads=1 still splits into multiple files, 
#which runs relatively quickly (~10min?)
THREADS=1

#then perform basecalling 
echo "alignments: ${ALIGNMENTS}"
echo "genome: ${GENOME}"
if [ ${PARALLEL} -ne 0 ]; then
	echo "parallel threads: ${THREADS}"
	echo "window size: ${PARALLEL}"
fi

if [ $PARALLEL -ne 0 ]; then
	if [ -d ${_TMPDIR} ]; then
			rm -rf ${_TMPDIR}
	fi

	mkdir ${_TMPDIR}
	bedtools makewindows -g ${GENOME}.fai -w ${PARALLEL} | bioawk -t '{print $1 ":" $2+1 "-" $3}' > ${_TMPDIR}/tmp.windows

	if [ $GVCF -eq 0 ]; then	
		parallel --jobs ${THREADS} gatk HaplotypeCaller -R ${GENOME} -I ${ALIGNMENTS} -O ${_TMPDIR}/tmp.{#}.vcf.gz --intervals {1} :::: ${_TMPDIR}/tmp.windows
		#then need to merge region'd vcfs
		bcftools concat -a --threads ${THREADS} -Oz ${_TMPDIR}/tmp.*.vcf.gz | bcftools sort -Oz -o ${_NAME}.snp.vcf.gz
	else
		parallel --jobs ${THREADS} gatk HaplotypeCaller -R ${GENOME} -I ${ALIGNMENTS} -O ${_TMPDIR}/tmp.{#}.vcf.gz -ERC GVCF --intervals {1} :::: ${_TMPDIR}/tmp.windows
		#then need to merge region'd vcfs
		bcftools concat -a --threads ${THREADS} -Oz ${_TMPDIR}/tmp.*.vcf.gz | bcftools sort -Oz -o ${_NAME}.g.vcf.gz
	fi
	rm -rf ${_TMPDIR}
else 
	if [ $GVCF -eq 0 ]; then
		gatk HaplotypeCaller -I ${ALIGNMENTS} -O ${_NAME}.snp.vcf.gz -R ${GENOME}
	else
		gatk HaplotypeCaller -I ${ALIGNMENTS} -O ${_NAME}.g.vcf.gz -R ${GENOME} -ERC GVCF
	fi
fi

###index VCF w/ GATK IndexFeatureFile for gVCF
if [ $GVCF -eq 1 ]; then 
	gatk IndexFeatureFile -I ${_NAME}*vcf.gz
fi

#after basecalling soft-filter variants for low quality (GQ>=20) and
#heterozygous genotypes
#only do this if not outputting gVCF
#if [ $GVCF -eq 0 ]; then
#	bcftools filter --threads ${THREADS} -s lowQual -i 'GQ>=10' ${_NAME}.snp.vcf.gz | bcftools filter --threads ${THREADS} -s Het -i 'GT="hom"' -o ${_NAME}.snp.soft-filter.vcf.gz -Oz
#	bcftools index ${_NAME}.snp.soft-filter.vcf.gz
#fi
