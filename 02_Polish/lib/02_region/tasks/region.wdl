version 1.0

task predata{    
    input{
        File hifi_fofn
        String workdir
        String sample
    }
    command <<<
        mkdir -p ${workdir}
        cd ${workdir}
        cat ~{hifi_fofn}|xargs >  ~{workdir}/new_hifidata
    >>>
    output{
        String  new_hifidata = workdir + '/new_hifidata'
    }
}

task minimap2_pb {
    input {
        File hifi_fofn
        String sample
        File software_json
        String workdir
        File genome

        Int cpu
        String memory
    }
    
    Array[Array[String]] hifidata = read_tsv(hifi_fofn)
    Object software = read_json(software_json)
    String minimap2 = software.minimap2
    String samtools = software.samtools

    command{
        set -ex
        cd ${workdir}
        if [-f '02.hifi_map_to_genome.sort.bam.done' ];then
            exit 0
        fi
        echo "
        set -ex
        cd ${workdir}
        ${minimap2} -ax map-pb -K 200M  -t ${cpu}  ${genome}  ${sep=" "  hifidata[0]} |${samtools} view  --threads 10 -T ${genome} -bS |${samtools} sort -t 8 -m 10g  -o  ${workdir}/02.hifi_map_to_genome.sort.bam
        ${samtools} index -@  5  ${workdir}/02.hifi_map_to_genome.sort.bam

        touch 02.hifi_map_to_genome.sort.bam.done

        " > minimap2_pb.sh
        bash  minimap2_pb.sh
    }
    output{
        File hifi_map_to_genome = workdir + '/02.hifi_map_to_genome.sort.bam'
    }
    runtime {
        cpu : cpu
        memory : memory
    }
}

task region {
    input {
        File minimap2_pb_bam
        File software_json
        String workdir
        File genome

        Int cpu
        String memory
    }
    
    Object software = read_json(software_json)
    String python3 = software.python3
    String script = software.script  

    command {
        set -ex        
        cd ${workdir}
        if [ -f '02.hifi_map_to_genome.sort.bam.nomap.done' ];then
            exit 0
        fi
        echo "
        set -ex
        cd ${workdir}
        ${python3}  ${script}/bam_region_depth.py  ${workdir}/02.hifi_map_to_genome.sort.bam  ${genome}  ${workdir}/02.hifi_map_to_genome.sort.bam  

        touch 02.hifi_map_to_genome.sort.bam.nomap.done
        " > region.sh
        bash region.sh
    }
    output {
        File nomap_fa = workdir + '/02.hifi_map_to_genome.sort.bam.nomap.false.fa'
        File gap_genome = workdir +'/02.hifi_map_to_genome.sort.bam.nomap.10000.fa'
    }
    runtime {
        cpu : cpu
        memory : memory
    }

}

