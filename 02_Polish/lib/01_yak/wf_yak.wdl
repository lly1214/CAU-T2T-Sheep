version 1.0

import "tasks/yak.wdl" as yak

workflow Yak {
    input {
        File software_json  #软件路径json格式
        String workdir  #工作目录
        String sample  #样本名称
        File sgs_fofn #二代数据绝对路径

        Int cpu = 10
        String memory = '50G'
    }
    call yak.shortyak_task as shortyak {
        input :
        software_json = software_json,
        workdir = workdir,
        sample = sample,
        sgs_fofn = sgs_fofn,
        cpu = cpu,
        memory = memory
    }
    call yak.longyak_task as longyak {
        input :
        software_json = software_json,
        workdir = workdir,
        sample = sample,
        sgs_fofn = sgs_fofn,
        cpu = cpu,
        memory = memory 
    }
}