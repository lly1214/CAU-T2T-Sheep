version 1.0

task nextpolish2_fit {
    input {
        File genome
        File software_json
        File workdir

        String cpu 
        String memory
    }

    Object software = read_json(software_json)
    String script = software.script

    command {
        set -ex
        cd ${workdir}
        if [ -f '05.asm.gap.closed.fa.hifipolish.filt.fbed.done' ];then
            exit 0
        fi
        echo "
        set -ex
        cd ${workdir}
        ${script}/nextPolish2 --out ${workdir}/05.asm.gap.closed.fa.hifipolish.filt.fbed  --out_pos  --thread 20 ${workdir}/02.hifi_map_to_genome.sort.bam   ${genome} ${workdir}/01.sr.short.yak  ${workdir}/01.sr.long.yak
        
        " > nextpolish2_fit.sh
        bash nextpolish2_fit.sh

    }
    output {
        File filtbed = workdir + '/05.asm.gap.closed.fa.hifipolish.filt.fbed'
        String wd4 = workdir 
    }
    runtime {
        cpu : cpu
        memory : memory
    }

}