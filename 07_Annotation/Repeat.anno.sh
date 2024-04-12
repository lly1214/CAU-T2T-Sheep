#!/bin/sh
#GMATA
perl gmata.pl -c default_cfg.txt -i hy.finallypolish.fasta

#TRF
TRF.py --TargetGenome hy.finallypolish.fasta --key hy --OutDir 02_TRF --rerun 2 Soft.conf

#MITE
MITE.py --TargetGenome hy.finallypolish.fasta --group 20 --percentage 0.1 --key hy --OutDir 03a_MITE --cpus 25 Soft.conf

#RepeatModeler
BuildDatabase -name hy -engine wublast hy.finallypolish.fasta.masked
RepeatModeler -pa 30 -database hy -engine wublast

#RepeatMasker
RepeatMasker.py --TargetGenome hy.finallypolish.fasta --key hy --OutDir FinalMask --cutf 15 --lib all.lib --species animal --cpus 40 --ssrgff 01_GMATA/hy.finallypolish.fasta.ssr.gff3 --TRFgff hy.TRF.gff3 Soft.conf