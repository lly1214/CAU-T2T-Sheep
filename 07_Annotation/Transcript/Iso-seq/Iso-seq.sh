#ccs calling
ccs --min-passes 1 --min-rq 0.9 --max-length 30000 --min-length 100 -j 10 split.5p--bc1002_3p.subreads.bam split.5p--bc1002_3p.ccs.bam

#Lima
lima split.5p--bc1002_3p.ccs.bam barcode.fa split.5p--bc1002_3p.fl.bam --isoseq --peek-guess -j 10

#Refine
isoseq3 refine --require-polya  split.5p--bc1002_3p.fl.bc1002_5p--bc1002_3p.bam  barcode.fa split.5p--bc1002_3p.flnc.bam

#Cluster
isoseq3 cluster flnc.fofn polished.bam --verbose --use-qvs -j 10 --singletons

#Minimap2
minimap2 --split-prefix split -ax splice -t 10 -uf --secondary=no -C5 hy.finallypolish.fasta polished.hq.fasta > minimap.sam
samtools sort -o minimap.sorted.sam minimap.sam

#Collapse
python rm_fusion_sam.py fusion.group.txt minimap.sorted.sam rm_fusion.sorted.sam
python collapse_isoforms_by_sam.py --input polished.hq.fasta -s rm_fusion.mapped.sorted.sam -o all
all.collapsed.rep.fa |awk -F '|' '{print $1}' > ISOSEQ_transcipt.fasta