# T2T-sheep1.0
We assembled a telomere-to-telomere (T2T) gap-free genome (2.85 Gb) from a ram (T2T-sheep1.0), which included a Y chromosome of 26.59 Mb, using ultralong ONT and PacBio HiFi reads together with Hi-C and Bionano data.  

The scripts for T2T-sheep1.0 are deposited in this website, including genome assembly, annotation, methylation, ChIP-seq, genome assessment, SD calling, centromere repeat identification, etc.  

The T2T genome polish pipeline was developed combining PacBio, ONT long reads, and short reads, and it employs NextPolish, NextPolish2, and other dependencies to improve the genome quality greatly.
## Descriptions 
### 01_Assembly
The 01_Assembly directory contains codes for trio binning, haploid genome assembly. The commands for Hi-C anchoring and assemblies of telomeres are also shwon in the directory.
### 02_Polish
The polish pipeline was written based on WDL. It involves 8 steps of processing, including (1) k-mer generation for the following NextPolish softwares, (2) determination of low-quality regions (LQRs, mostly referring to the gaps) and high-quality regions (HQRs), (3) LQRs mapped for their chromosomal coordinates and orientations, (4) polishing LQRs by NextPolish, (5) polishing LQRs by NextPolish2, (6) polishing high coverage regions (HQRs) by NextPolish2, (7) merging LQRs and HQRs and (8) final polish for the whole genome by NextPolish2.
#### Dependencies
* Python 3.6+
* java 11.0+
* Cromwell  
Other dependent softwares are included in the \bin\bin directory
#### Input
Low-quality genome fasta file to be polished is required.
The raw data (fastq) includes all the short and long reads, such as, PacBio long reads, ONT long reads, and short reads (MGI or Illumina).
#### Usage
1. Provide the environment paths of the required softwares in the software.json file.
``` 
{
	"yak":"/02_Polish/bin/bin/yak",
	"minimap2":"/02_polish/bin/bin/minimap2",
	"samtools":"/02_Polish/bin/bin/samtools",
	"python3":"/Path/Python-3.6.3/bin/python3",
	"script":"/02_Polish/bin/bin/",
	"Workflow":"/02_Polish/bin/"
}  
```
2. The running information is provided in the run.json file, including HiFi reads (hifi.fofn), memeory in the workstation, ONT reads (lgs.fofn), genome file to be polished (test.fa), sample name (test), short reads (sgs.fofn), cpus in the workstation (cpu).
```
{
  "nextPolish2.hifi_fofn": "/PATH/hifi.fofn",
  "nextPolish2.workdir": "/PATH/test/",
  "nextPolish2.memory": "60G",
  "nextPolish2.lgs_fofn": "/PATH/lgs.fofn",
  "nextPolish2.genome": "/Path/test.fa",
  "nextPolish2.sample": "test",
  "nextPolish2.software_json": "/PATH/software.json",
  "nextPolish2.sgs_fofn": "/PATH/sgs.fofn",
  "nextPolish2.cpu": 20
}
```
3. Run the polish workflow in the terminal.
```
bash run.sh
```   
#### Output
The polished genome .fa file is finally generated, with the consensus quality value (QV) of the polished genome is also provided.
  
### 03_Assessment
The 03_Assessment directory contains codes for colinearity analysis compared with current reference, BUSCO evaluation, Quality Value (QV) assessment and coverage calculating of raw data (HiFi, ONT, NGS).

### 04_SD
The 04_SD directory contains codes for SD calling and circos plot.

### 05_Methylation
The 05_Methylation directory contains codes for methylation calling based on HiFi and ONT sequencing data. 

### 06_ChIP-seq 
The 06_ChIP-seq directory contains codes for mapping and processing ChIP-seq data.

### 07_Annotation
The 07_Annotation directory contains codes for repeat annotation and gene annotation based on the three strategies of transcript-based, *de novo*, and homolog-based predictions.

### 08_Centromeres
The 08_Centromeres directory contains codes for (1) centromeric sequence similarity calculation and heatmap visualization using [StainedGlass](https://github.com/mrvollger/StainedGlass), (2) calculation of the sequence complexity and entropy using [NeSSie](https://github.com/projectnessie/nessie) and (3) identification of satellite repeats using [SRF](https://github.com/lh3/srf).

### 09_OrthoFinder 
The 09_OrthoFinder directory contains codes for identifying orthologous gene groups using [OrthoFinder](https://github.com/davidemms/OrthoFinder).

