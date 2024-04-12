version 1.0

task shortyak_task{ # calculate 21kmer
    input {
        File software_json
        String workdir
        String sample 
        File sgs_fofn

        String cpu
        String memory
    }

    Object software = read_json(software_json)
    String yak = software.yak
    Array[Array[String]] sgs_data = read_tsv(sgs_fofn)

    command {
        set -ex
        cd ${workdir}
        if [ -f '01.shortyak.done' ]; then
            exit 0
        fi
        echo " 
        cd ${workdir}
        ${yak}  count  -k 21 -t ${cpu} -o ${workdir}/01.sr.short.yak  <(zcat ${sep=" " sgs_data[0]})   <(zcat ${sep=" " sgs_data[1]})
        
        if [ $? -eq 0 ];then
            touch ${workdir}/01.shortyak.done
        fi
            touch ${workdir}/shortyak.err            
        " > shortyak.sh
        bash shortyak.sh      
    }
    output {
        File shortyak = workdir + '/01.sr.short.yak'
    }
    runtime {
        cpu : cpu
        memory : memory
    }
}

task longyak_task{ # calculate 31kmer
    input {
        File software_json
        String workdir
        String sample 
        File sgs_fofn

        String cpu
        String memory
    }   

    Object software = read_json(software_json)
    String yak = software.yak
    Array[Array[String]] sgs_data = read_tsv(sgs_fofn)

    command {
         if [ -f '01.longyak.done' ];then
            exit 0
        fi
        echo " 
        cd ${workdir}
        ${yak}  count  -k 31 -t ${cpu} -o ${workdir}/01.sr.long.yak    <(zcat ${sep=" " sgs_data[0]})   <(zcat ${sep=" " sgs_data[1]})
        if [ $? -eq 0 ];then
            touch ${workdir}/01.longyak.done
        fi
            touch ${workdir}/longyak.err    
        " > longyak.sh
        bash longyak.sh        
    }
    output {
        File longyak = workdir + '/01.sr.long.yak'
        String wd1 = workdir
    }
    runtime {
        cpu : cpu
        memory : memory 
    }
}


