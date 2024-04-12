version 1.0

task nextpolish2_nofit {
    input {
        File software_json
        File nomap_fa
        File workdir

        String cpu 
        String memory
    }

    Object software = read_json(software_json)
    String script = software.script

    command {
        set -ex
        cd ${workdir}
        if [ -f '06.asm.gap.closed.fa.hifipolish.nofilt.fbed.done' ];then
            exit 0
        fi
        echo "
        set -ex
        cd ${workdir}
        ${script}/nextPolish2 --out ${workdir}/06.asm.gap.closed.fa.hifipolish.nofilt.fbed  --out_pos  --thread 20 --use_secondary  --min_map_qual -1  ${workdir}/02.hifi_map_to_genome.sort.bam    ${nomap_fa}  ${workdir}/01.sr.short.yak  ${workdir}/01.sr.long.yak
        
        " > nextpolish2_nofit.sh
        bash nextpolish2_nofit.sh

    }
    output {
        File nofiltbed = workdir + '/06.asm.gap.closed.fa.hifipolish.nofilt.fbed'
        String wd5 = workdir
    }
    runtime {
        cpu : cpu
        memory : memory
    }

}