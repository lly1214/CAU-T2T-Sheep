#!/bin/sh
name=$1
sequence=$2
database=mammalia_odb10
mode=geno  

python busco --cpu 16 --mode $mode --force --lineage_dataset $database --offline --in $sequence --out $name.out  --out_path $name.work.out