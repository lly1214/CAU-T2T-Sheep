version 1.0

task merge {
    input {
        File software_json
        File genome
        File sgs_fofn
        File hifi_fofn
        File lgs_fofn
        String workdir

        String cpu
        String memory
    }

    Object software = read_json(software_json)
    String script = software.script
    String python3 = software.python3
    String Workflow = software.Workflow

    command {
        set -ex
        cd ${workdir}
        if [ -f '07.asm.gap.closed.fa.polish.merge.fasta.done' ];then
            exit 0
        fi
        echo "
        set -ex
        cd ${workdir}

        ${python3} ${Workflow}/nextp2w.py  --sr ${sgs_fofn}  --ont ${lgs_fofn} --hifi ${hifi_fofn} --genome ${genome} --workdir ${workdir} --thread ${cpu}

        " > merge.sh
        bash merge.sh
    }

    output {
        File mergefasta = workdir + '/07.asm.gap.closed.fa.polish.merge.fasta'
    }

    runtime {
        cpu : cpu
        memory : memory
    }

}
