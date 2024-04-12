samtools sort -n NGS/STAR/Aligned.out.bam -@ 5 -o Aligned.out.ssn.bam
Augustus/latest/bin/filterBam --uniq --in Aligned.out.ssn.bam --out Aligned.out.ssf.bam
samtools sort Aligned.out.ssf.bam -@ 5 -o Aligned.out.ssfs.bam
Augustus/latest/bin/bam2hints --intronsonly --in=Aligned.out.ssfs.bam --out=introns.gff
BRAKER/scripts/filterIntronsFindStrand.pl hy.finallypolish.fasta introns.gff --score > hint.gff