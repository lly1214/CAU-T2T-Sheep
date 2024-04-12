#!/bin/sh
#21-mer database building
for i in r1 r2 ;
do
    meryl k=21  threads=15   memory=100g  count output read${i}.meryl HU-3095.clean.${i}.fastq.gz
done

meryl union-sum output   HU.k21.meryl  read*.meryl

#quality value assessment
merqury.sh   HU.k21.meryl   hy.finallypolish.fasta   hy