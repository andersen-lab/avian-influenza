#!/bin/bash

read1=$1;
read2=$2;
n=$(basename $read1);
n=${n/_1.fastq/};

index=/alfheim/h5n1/ref/EPI_ISL_19032063

echo $read1 $read2 ${n}

# Align
~/code/bwa/bwa mem -t 5 $index $read1 $read2 | samtools view -F 4 -b | samtools sort -o ${n}.bam 


