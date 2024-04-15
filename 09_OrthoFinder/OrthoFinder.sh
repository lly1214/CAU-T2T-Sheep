#!/bin/sh
orthofinder -f  hy_pep     -S  blast  -M msa -A mafft  -T  raxml  -t 10  -a 10 

#You should put pep.fasta files of all the related species to a directory named "hy_pep"