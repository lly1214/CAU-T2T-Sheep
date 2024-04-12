version 1.0

task finally_polish {
    input {
        File software_json
        File hifi_fofn
        File mergefasta
        String workdir 

        String cpu
        String memory
    }

    Array[Array[String]] hifidata = read_tsv(hifi_fofn)
    Object software = read_json(software_json)
    String script = software.script
    String python3 = software.python3
    String samtools = software.samtools
    String minimap2 = software.minimap2

    command {
        set -ex
        cd ${workdir}
        if [ -f '08.merge.fasta.sort.bam.done && 08.finally_polish.done' ] 
            exit 0
        fi
        echo "
        cd ${workdir}
        #minimap2 hifi 
        ${minimap2} -ax map-hifi  -t ${cpu}  ${mergefasta}  ${sep=" "  hifidata[0]} |${samtools} view  --threads 10 -T ${mergefasta} -bS |${samtools} sort -t 8 -m 10g  -o  ${workdir}/08.merge.fasta.sort.bam
        ${samtools} index -@ 10  08.merge.fasta.sort.bam
        touch 08.merge.fasta.sort.bam.done 

        #finally_polish 
        ${script}/nextPolish2  --out  ${workdir}/08.finally_polish.fasta  --thread ${cpu}  ${workdir}/08.merge.fasta.sort.bam  ${mergefasta}  ${workdir}/01.sr.long.yak  ${workdir}/01.sr.short.yak
        touch 08.finally_polish.done

        " > finally_polish.sh
        bash finally_polish.sh
    }

    output {

    }
    runtime {
        cpu : cpu
        memory : memory
    }

}