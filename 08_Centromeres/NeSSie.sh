name='HuYang'

for i in ${workdir}/*fasta
do

idx=`basename $i .fasta`
nessie -I $i -O ${name}_${idx}_nessie1000_8.out -E -l 1000  -s 8

done

#plot
python SDT_pdf_ForNessie.py -i $nessiseresult  -Low 0.50  -Xmax ${XmaxNum}  -Xparts  10  -Pl 30 -Xid  ${idx}   -size 0.00001   -o ./${idx}_plot/\${idx}_0.50_1_nessie1000_8.png