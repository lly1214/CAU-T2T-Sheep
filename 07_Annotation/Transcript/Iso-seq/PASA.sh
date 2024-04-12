#get transcript
gffread -w ngs_transcript.fasta -g hy.finallypolish.fasta NGS/StringTie/merged.gtf
ln -s Collapse/ISOSEQ_transcipt.fasta tgs_transcript.fasta

#PASA
python3 PASA.py --genome hy.finallypolish.fasta --transngs ngs_transcript.fasta --transtgs tgs_transcript.fasta --key hy --cores 5 --workdir PASA