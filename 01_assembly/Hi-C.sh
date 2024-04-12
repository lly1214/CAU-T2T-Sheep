#!/bin/sh
#quality control
fastp -i HU-3095.R1.fq.gz -I HU-3095.R2.fq.gz -o HU-3095.clean_R1.fq.gz -O HU-3095.clean_R2.fq.gz -R HU-3095 -h HU-3095.QC.html -j HU-3095.QC.json -s 40 -w 4

#mapping
bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder --phred33-quals -p 4  --no-sq --un trim/0005.HU-3095.R2.unmap.fastq --rg-id BMG --rg SM:0005.HU-3095.R2 -x /hi-c/midfile/HU-3095.Contig.fasta -U /hi-c/01_QC/0005.HU-3095.clean_R2.fq.gz | samtools view -F 4 -@ 2 -t /hi-c/midfile/HU-3095.Contig.fasta.fai -bS -> tmp/0005.HU-3095.R2.bam

#valid pairs generation (HiC-pro)
mergeSAM.py -q 0 -t -v -f 0007.HU-3095.R1.Merge.sort.bam -r 0007.HU-3095.R2.Merge.sort.bam -o 0007.HU-3095.paired.sort.bam
mapped_2hic_fragments.py -f HU-3095.Contig.dpnii.bed -r 0001.HU-3095.paired.sort.bam -v -S -t 100 -m 100000000 -s 100 -l 700 -a -o .

#anchor
Lachesis HU-3095.ini 
