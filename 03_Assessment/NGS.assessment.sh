#!/bin/sh
#quality control
fastp -i HU-3095_raw_1.fq.gz -I HU-3095_raw_2.fq.gz \
-o HU-3095.clean.r1.fastq -O HU-3095.clean.r2.fastq \
--thread 15 -n 0 -f 5 -F 5 -t 5 -T 5 \

#mapping
bwa mem -t 20  -M  -k 32   hy.finallypolish.fasta  HU-3095.clean.r1.fastq.gz  HU-3095.clean.r2.fastq.gz   | samtools view --threads 5 -bS | samtools sort --threads 5  -m  10G -o  hy.sort.bam

#coverage calculating
bedtools makewindows -g chr.len -w 200000 -s 200000 >ref.bed    
awk '{s++;print$1"\t"$2+1"\t"$3}' ref.bed >ref.new.bed 
bamdst -p ref.new.bed -o  ${outdir}   ${bam}

#bam_stat
samtools  stats -d  --threads 10  hy.sort.bam > hy.ngs.stat
samtools  flagstat --threads  10  hy.sort.bam > hy.ngs.flagstat