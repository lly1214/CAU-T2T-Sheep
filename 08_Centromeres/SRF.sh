#!bin/bash

mkdir -p tmp_dir

# set -ci to approximately 10-fold over the average coverage
kmc -fa -k151 -t16 -ci100 -cs100000  @hifi.fofn  count.kmc tmp_dir
kmc_dump  count.kmc count.txt

# assemble satellite DNA
srf -p prefix count.txt > srf.fa


# analyze
minimap2 -c -N 1000000 -f 1000 -r 100,100 -t 20  <(minimap2/k8-Linux  srfutils.js enlong  srf.fa )  hy.finallypolish.fasta   > srf-aln.paf
minimap2/k8-Linux  srfutils.js paf2bed  srf-aln.paf > srf-aln.bed   # filter and extract non-overlapping regions
minimap2/k8-Linux	 srfutils.js bed2abun srf-aln.bed > srf-aln.len  # calculate abundance of each srf contig
