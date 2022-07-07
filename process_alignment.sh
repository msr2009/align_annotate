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
	echo "syntax: process_alignment.sh --input BAM/SAM --threads THREADS"
	echo 
}


#check for appropriate argument
#if not there, kill and print help
#if [ $# -eq 1 ]
#then
#	if [[ ($1 == *.sam) || ($1 == *.bam) ]]
#	then
#		ALIGNMENTS="$1"
#	else
#		echo "ERROR: incorrect input format. requires .sam or .bam"
#		echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#		echo
#		HELP
#		exit
#	fi
#else		
#	HELP
#	exit
#fi

THREADS=1

if [ $# -eq 0 ]
then
	HELP
	exit
fi	

while [ $# -gt 0 ]; do
    case "$1" in
		-i|--input)
			ALIGNMENTS="$2"
			shift 2
			;;
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-h|--help) #HELP ME PLEASE!
			HELP
			exit;;
		*)
			HELP
            exit;;
    esac
done

FEWER_THREADS=$(( $THREADS - 1 ))

_NAME="${ALIGNMENTS%%.*}"

#if the file isn't a bam already, then make one
if [[ ${ALIGNMENTS} == *.sam ]] 
then
	samtools view -bS ${_NAME}.bam
fi

#then do all the processing	
samtools collate -@ ${FEWER_THREADS} -O ${_NAME}.bam | samtools fixmate -@ ${FEWER_THREADS} -m - - | samtools sort -@ ${THREADS} - | samtools markdup -@ ${FEWER_THREADS} -s - ${_NAME}.srt.rmdup.bam
samtools index -@ ${FEWER_THREADS} ${_NAME}.srt.rmdup.bam ${_NAME}.srt.rmdup.bam.bai
