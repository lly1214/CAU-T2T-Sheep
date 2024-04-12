version 1.0

task nextPolish1 {
    input {
        File software_json 
        String workdir
        String sample 
        File sgs_fofn
        File lgs_fofn
        File hifi_fofn
        File genome

        String cpu
        String memory 
    }
    Array[Array[String]] ontdata = read_tsv(lgs_fofn)
    Array[Array[String]] sgs_data = read_tsv(sgs_fofn)
    Object software = read_json(software_json)
    String minimap2 = software.minimap2
    String samtools = software.samtools
    String python3 = software.python3
    String script = software.script
    String Workflow = software.Workflow

    command {
        set -ex
        cd ${workdir}
        if [ -f '04.hifi_map_to_genome.sort.bam.ont.ref.fa.nextpolish.fa.done' ];then
            exit 0
        fi
        echo " 
        #nextPolish1 step1
        set -ex
        #get gap_genome
        cd ${workdir}
        #${python3} ${script}/generate_ont_gap_ref.py  --hifi_map_ref 03.hifi_map_to_gap.sort.bam  --hifi_gap_filt_ref  02.hifi_map_to_genome.sort.bam.nomap.10000.fa --ont_gap_ref 02.hifi_map_to_genome.sort.bam.nomap.false.fa --out_file 04.hifi_map_to_genome.sort.bam.ont.ref.fa
        ${python3} ${Workflow}/nextp2w.py  --sr ${sgs_fofn}  --ont ${lgs_fofn} --hifi ${hifi_fofn} --genome ${genome} --workdir ${workdir} --thread ${cpu}

        touch  04.hifi_map_to_genome.sort.bam.ont.ref.fa.done
        "> nextPolish1.sh 
        bash nextPolish1.sh
    }

    output {
        File nextPolish_genome  = workdir + '/04.hifi_map_to_genome.sort.bam.ont.ref.fa.nextpolish.fa'
        String wd3 = workdir
    }
    
    runtime {
        cpu : cpu
        memory : memory

    }
}
