#!/bin/sh
#mapping
nanopolish index -d 20220706-UNL310-P6-PAK30511-fast5/ 20220706-UNL310-P6-PAK30511-sup.pass.fastq.gz 
minimap2  -K 200M   --secondary=no -t 20 -ax  map-ont   hy.finallypolish.fasta   20220706-UNL310-P6-PAK30511-sup.pass.fastq.gz  | samtools view --threads 10 -T  hy.finallypolish.fasta   -bS | samtools sort --threads 8 -m 1G -o  20220706-UNL310-P6-PAK30511-sup.sort.bam 
samtools index  -@ 10  20220706-UNL310-P6-PAK30511-sup.sort.bam

#coverage calculating
bedtools makewindows -g chr.len -w 200000 -s 200000 >ref.bed    
awk '{s++;print$1"\t"$2+1"\t"$3}' ref.bed >ref.new.bed 
bamdst -p ref.new.bed -o  ${outdir}   ${bam}

#bam_stat
samtools  stats -d  --threads 10  hy.sorted.bam  > hy.ont.stat
samtools  flagstat --threads  10  hy.sorted.bam  > hy.ont.flagstat