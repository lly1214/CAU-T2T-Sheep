version 1.0

task Predata{
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

task map2gap {
    input {
        File hifi_fofn
        String sample
        File software_json
        String workdir
        String gap_genome

        Int cpu
        String memory
    }
    
    Array[Array[String]] hifidata = read_tsv(hifi_fofn)
    Object software = read_json(software_json)
    String minimap2 = software.minimap2
    String samtools = software.samtools
    String python3 = software.python3
    String script = software.script

    command{
        set -ex
        cd ${workdir}
        if [-f '03.hifi_map_to_gap.sort.bam.done'&&'03.hifi_map_to_gap.sort.bam.filt.done' ];then
            exit 0
        fi
        echo "
        set -ex
        #map2gap genome
        cd ${workdir}
        ${minimap2} -ax map-pb  -t ${cpu}  ${gap_genome}  ${sep=" "  hifidata[0]} |${samtools} view  --threads 10 -T ${gap_genome} -bS |${samtools} sort -t 8 -m 10g  -o  ${workdir}/03.hifi_map_to_gap.sort.bam
        ${samtools} index -@  5  ${workdir}/03.hifi_map_to_gap.sort.bam
        touch 03.hifi_map_to_gap.sort.bam.done

        #filterbam
        ${python3} ${script}/bam_filter.py --min_map_depth 3 --genome  ${gap_genome}   ${workdir}/03.hifi_map_to_gap.sort.bam   ${workdir}/03.hifi_map_to_gap.sort.bam.filt
        ${samtools} index ${workdir}/03.hifi_map_to_gap.sort.bam.filt.bam        
        touch  03.hifi_map_to_gap.sort.bam.filt.done
        " > minimap2_gap.sh
        bash  minimap2_gap.sh
    } 
    output{
        File filt_bam = workdir + '/03.hifi_map_to_gap.sort.bam.filt.bam'
        File filt_fa = workdir + '/03.hifi_map_to_gap.sort.bam.filt.fa'
    }
    runtime {
        cpu : cpu
        memory : memory
    }
}

task polish_gap{
    input {
        String workdir
        File filt_fa
        File filt_bam
        File software_json

        Int cpu
        String memory
    }

    Object software = read_json(software_json)
    String minimap2 = software.minimap2
    String samtools = software.samtools
    String python3 = software.python3
    String script = software.script

    command{
        set -ex
        cd ${workdir}
        if [ -f '03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa.fbed.done' ];then
            exit 0
        fi
        echo "
        set -ex
        cd ${workdir}
        #nextPolish2 round1
        rm -r -f  ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa
        ${script}/nextPolish2  -t ${cpu}  -o ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa    ${workdir}/03.hifi_map_to_gap.sort.bam.filt.bam  ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa  ${workdir}/01.sr.short.yak  ${workdir}/01.sr.long.yak   
        #
        ${minimap2} -ax asm20 -Y -t ${cpu}  ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa  ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa  |${samtools} sort -t 20 -o  ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa.bam
        #
        ${python3} ${script}/bam2fbed.py ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa.bam > ${workdir}/03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa.fbed

        touch 03.hifi_map_to_gap.sort.bam.filt.fa.hifipolish.filt.fa.fbed.done
        " > polish_gap.sh
        bash polish_gap.sh
    }
    output {
        String wd2 = workdir
    }
    runtime {
        cpu : cpu
        memory : memory
    }
}