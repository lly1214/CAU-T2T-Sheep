#!bin/bash

#mapping
minimap2   -t 10  -cx asm10  hy.finalpolish.fa Ramb_v3.0.fa >  hy-ram3.paf

#plot
paf2dotplot.r   -b  hy-ram3.paf