#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#REQUIRES: samtools

HELP(){
	echo "process_alignment.sh"
	echo "wrapper script for samtools commands to process alignments"
	echo "takes positional argument of BAM or SAM file containing alignments"
	echo
	echo "syntax: process_alignment.sh [BAM/SAM]"
}


#check for appropriate argument
#if not there, kill and print help
if [ $# -eq 1 ]
then
	if [[ ($1 == *.sam) || ($1 == *.bam) ]]
	then
		ALIGNMENTS="$1"
	else
		echo "ERROR: incorrect input format. requires .sam or .bam"
		echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		echo
		HELP
		exit
	fi
else		
	HELP
	exit
fi

_NAME="${ALIGNMENTS%%.*}"

#if the file isn't a bam already, then make one
if [[ ${ALIGNMENTS} == *.sam ]] 
then
	samtools view -bS ${_NAME}.bam
fi

#then do all the processing	
samtools collate -O ${_NAME}.bam | samtools fixmate -m - - | samtools sort - | samtools markdup -s - ${_NAME}.srt.rmdup.bam
samtools index ${_NAME}.srt.rmdup.bam ${_NAME}.srt.rmdup.bam.bai
