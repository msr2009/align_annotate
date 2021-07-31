#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#REQUIRES: bcftools

HELP(){
        echo "call_variants.sh"
        echo "wrapper script for bcftools mpileup and call commands"
        echo "takes two positional arguments:"
		echo "	1) BAM or SAM file containing alignments"
		echo "  2) FASTA file for reference genome"
		echo 
        echo "syntax: call_variants.sh [BAM/SAM] [FASTA]"
}

###this seems excessive but it works? 

if [ $# -eq 2 ] 
then 
		while [ $# -gt 0 ]; do
			case "$1" in 
					*.sam|*.bam)
							ALIGNMENTS="$1"
							shift 1
							;;
					*.fa|*.fasta|*.fsa)
							GENOME="$1"
							shift 1
							;;
			esac
		done
else
		HELP
		exit
fi

_NAME="${ALIGNMENTS%%.*}"

#then perform basecalling 
bcftools mpileup -Ou -f ${GENOME} ${ALIGNMENTS} | bcftools call -mv -Ov -o ${_NAME}.vcf
