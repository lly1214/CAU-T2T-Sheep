#!/bin/sh
#mapping
bowtie2   -p 20   --very-sensitive --no-mixed --no-discordant -k 10   -x  hy.finallypolish.fasta    -1  all.CENP.R1.fq.gz  -2  all.CENP.R2.fq.gz  | samtools sort -O bam -@ 10 -o - >  hy.CENP.bam

#coverage
bedtools makewindows -g chr.len -w 10000 -s 10000 >ref.bed    
awk '{s++;print$1"\t"$2+1"\t"$3}' ref.bed >ref.new.bed 
bamdst -p ref.new.bed -o  ${outdir}   ${bam}