#!/bin/bash

insam=$1
inbam=$2
insort=$3
index=$4
contig=$5
cStart=$6
cEnd=$7
outbam=$8

cd /home/ubuntu/proc/sjspence/170105_PSE
#samtools view -b $insam > $inbam
#samtools sort -o $insort $inbam
#samtools index $insort > $index
samtools view $insort \"${contig}:${cStart}-${cEnd}\" > $outbam
