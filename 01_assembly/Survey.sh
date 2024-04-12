#!/bin/sh
#kmc 
/v3.1.1/kmc  -k21 -t16 -m64 -ci1 -cs10000000  @hifi.fofn  kmcdb tmp
#histo
kmc/v3.1.1/kmc_tools transform kmcdb histogram sample.histo -cx10000000 
#R
R/bin/Rscript /genomescope2.0/genomescope.R -i  sample.histo   -k 21  -o ./ -m -1

##write the path of your hifi data to hifi.fofn  