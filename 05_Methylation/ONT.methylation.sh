#!/bin/sh
#index
nanopolish index -d 20220706-UNL310-P6-PAK30511-fast5/ 20220706-UNL310-P6-PAK30511-sup.pass.fastq.gz

#mapping
minimap2  -K 200M   --secondary=no -t 20 -ax  map-ont   hy.finallypolish.fasta   20220706-UNL310-P6-PAK30511-sup.pass.fastq.gz  | samtools view --threads 10 -T  hy.finallypolish.fasta   -bS | samtools sort --threads 8 -m 1G -o  20220706-UNL310-P6-PAK30511-sup.sort.bam 
samtools index  -@ 10  20220706-UNL310-P6-PAK30511-sup.sort.bam

#methylation calling
nanopolish call-methylation -t 20 -r  20220706-UNL310-P6-PAK30511-sup.pass.fastq.gz  -b 20220706-UNL310-P6-PAK30511-sup.sort.bam   -g hy.finallypolish.fasta   > 20220706-UNL310-P6-PAK30511-sup.methylation_calls.tsv