#!/bin/bash

read1=$1;
read2=$2;
n=$(basename $read1);
n=${n/_1.fastq/};

index=/alfheim/h5n1/ref/EPI_ISL_19032063

echo $read1 $read2 ${n}


# Index
samtools index ${n}.bam

genes=(
"PB1|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"MP|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"NP|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"PB2|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"PA|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"HA|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"NS|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
"NA|A/dairy_cattle/Texas/24-008749-002-v/2024|EPI_ISL_19032063|2024-03-20|2.3.4.4b"
)

# Call consensus
for gene in "${genes[@]}"; do
	echo $gene
	prefixGene=$(echo $gene | cut -f 1 -d \|)
	echo $gene $prefixGene
	samtools mpileup -r $gene -A -d 0 -aa -Q 0 ${n}.bam | ~/code/ivar/build/bin/ivar consensus -p cns_seqs/${n}_${prefixGene}_cns -t 0.5 -m 1
done


