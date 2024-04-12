#!/bin/sh
#blast
makeblastdb  -in  fasta  -input_type fasta -dbtype nucl -title CL -out CL_db
blastn  -db  CL_db  -query  ref.fasta
#hifiasm assembly
hifiasm -t 20 -o telomere  all.hifi.fa
#ragtag
ragtag.py correct -t 10  --aligner  minimap2   ${reference}      ${query}
ragtag.py scaffold  -t 10  --aligner  minimap2    ${reference}   ./ragtag_output/ragtag.correct.fasta  -C
