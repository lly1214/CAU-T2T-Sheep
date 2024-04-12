#!/bin/sh
#mapping
samtools fasta  HU-3095.ccs.bam > HU-3095.ccs.fa 
pbmm2  align  --num-threads 24 --preset CCS  --sample hy  --log-level INFO --sort --unmapped   hy.finallypolish.fasta   input.fofn   hy.align.bam

#coverage calculating
bedtools makewindows -g chr.len -w 200000 -s 200000 >ref.bed    
awk '{s++;print$1"\t"$2+1"\t"$3}' ref.bed >ref.new.bed 
bamdst -p ref.new.bed -o  ${outdir}   ${bam}

#bam_stat
samtools  stats -d  --threads 10  hy.align.bam > hy.hifi.stat
samtools  flagstat --threads  10 hy.align.bam > hy.hifi.flagstat