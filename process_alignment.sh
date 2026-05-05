#!/bin/bash
set -o nounset
set -o errexit
#set -o pipefail

#REQUIRES: samtools

HELP(){
	echo "process_alignment.sh"
	echo "wrapper script for samtools commands to process alignments"
	echo "takes positional argument of BAM or SAM file containing alignments"
	echo
	echo "syntax: process_alignment.sh --input BAM/SAM --threads THREADS"
	echo 
}

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

case ${ALIGNMENTS} in *.sam)
		samtools view -bo ${_NAME}.bam ${_NAME}.sam
esac

#apparently I can't do this for posix-ness
#if [[ ${ALIGNMENTS} == *.sam ]] 
#then
#	samtools view -bS ${_NAME}.bam
#fi

#then do all the processing

#old processing 
#samtools sort -@ ${THREADS} ${_NAME}.bam | samtools markdup -@ ${FEWER_THREADS} -s - ${_NAME}.srt.rmdup.bam
#
#new processing w/ collate and fixmate should work:
#samtools collate -@ ${FEWER_THREADS} -Ou ${_NAME}.bam | samtools fixmate -@ ${FEWER_THREADS} -m - - | samtools sort -@ ${THREADS} - | samtools markdup -@ ${FEWER_THREADS} -s - ${_NAME}.srt.rmdup.bam
#
#but it doesn't, so we need to make some temporary files
#

echo `date +”%m/%d/%Y %H:%M”` collating reads
samtools collate -@ ${FEWER_THREADS} ${_NAME}.bam ${_NAME}.coll
echo `date +”%m/%d/%Y %H:%M”` fixing mates
samtools fixmate -@ ${FEWER_THREADS} -m ${_NAME}.coll.bam ${_NAME}.mate.bam
echo `date +”%m/%d/%Y %H:%M”` sorting, removing duplicates
samtools sort -@ ${THREADS} ${_NAME}.mate.bam | samtools markdup -@ ${FEWER_THREADS} -s - ${_NAME}.srt.rmdup.bam
echo `date +”%m/%d/%Y %H:%M”` indexing BAM
samtools index -@ ${FEWER_THREADS} ${_NAME}.srt.rmdup.bam ${_NAME}.srt.rmdup.bam.bai
chmod a+r ${_NAME}.srt.rmdup.bam*
