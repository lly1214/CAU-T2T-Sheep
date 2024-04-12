stringtie -o stringtie.gtf -p 8 Aligned.out.bam
echo "stringtie.gtf" > gtf.list
stringtie --merge -p 8 -o merged.gtf gtf.list
