#!/bin/sh
#quality control

fastp -i /HUF/HUF_raw_1.fq.gz -I /HUF/HUF_raw_2.fq.gz \
-o HUF_r1.fastq -O HUF_r2.fastq \
--thread 15 -n 0 -f 5 -F 5 -t 5 -T 5 \

fastp -i /HUM/HUM_raw_1.fq.gz -I /HUM/HUM_raw_2.fq.gz \
-o HUM_r1.fastq -O HUM_r2.fastq \
--thread 15 -n 0 -f 5 -F 5 -t 5 -T 5 \

#yak generation
yak count -t 20  -o  HUF.yak  <(zcat /HUF/HUF_r1.fastq.gz )  <(zcat /HUF/HUF_r2.fastq.gz )
yak count -t 20  -o  HUM.yak  <(zcat /HUM/HUM_r1.fastq.gz )  <(zcat /HUM/HUM_r2.fastq.gz )

#HiFiasm trio binning assembly
hifiasm  -o HU  -t 20   -1 HUF.yak -2 HUM.yak  HU3095.hifi.ccs.fa