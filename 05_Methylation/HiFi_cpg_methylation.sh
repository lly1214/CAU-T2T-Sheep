#!/bin/sh
#01_pbmm2
pbmm2  align  --num-threads 24 --preset CCS  --sample hy  --log-level INFO --sort --unmapped   hy.finallypolish.fasta   input.fofn   hy.align.bam

#02_pb-CpG-tools
pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores  --bam  hy.align.bam \
	 --ref   hy.finallypolish.fasta    --output-prefix  hy \
	 --model pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite   --threads  24   --min-mapq 10