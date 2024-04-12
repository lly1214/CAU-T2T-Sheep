version 1.0

import "./lib/01_yak/tasks/yak.wdl" as yak
import "./lib/02_region/tasks/region.wdl" as region
import "./lib/03_map_gap/tasks/map_gap.wdl" as map_gap
import "./lib/04_nextpolish1/tasks/task.wdl"  as nextpolish1
import "./lib/05_nextpolish2_fit/tasks/nextpolish2_fit.wdl" as nextpolish2fit
import "./lib/06_nextpolish2_nofit/tasks/nextpolish2_nofit.wdl"  as nextpolish2nofit
import "./lib/07_merge/tasks/mergefasta.wdl" as mergefasta
import "./lib/08_merge_polish/tasks/finallypolish.wdl"  as finallypolish
 
workflow nextPolish2 {

    input {
        File software_json  #JSON format for software path
        String workdir  #workdir
        String sample  #sample name
        File sgs_fofn #The absolute path of the NGS data
        File hifi_fofn #The absolute path of the HiFi data
        File lgs_fofn #The absolute path of the ONT data
        File genome #The path of genome .fa 

        Int cpu = 20
        String memory = '60G'
    }
   
    parameter_meta {
        workdir:"main workdir"
        software_json:"software path"
        sample:"sample name"
        sgs_fofn:"ngs data path"
        hifi_fofn:"hifi data path"
        lgs_fofn:"ont data path"
        genome:"genome path"
        cpu:"process number"
        memory:"mermory"
    }
    #steps
    #step1
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
    #step2
    call region.predata as predata {
            input :
            workdir = workdir,
            sample = sample,
            hifi_fofn = hifi_fofn
    }
    call region.minimap2_pb as minimap2_pb {
        input :
            hifi_fofn = predata.new_hifidata,
            sample = sample,
            software_json = software_json,
            workdir = longyak.wd1,
            genome = genome,
            cpu = cpu,
            memory = memory
    }
    call region.region as region {
        input :
            minimap2_pb_bam = minimap2_pb.hifi_map_to_genome,
            software_json = software_json,
            workdir = workdir,
            genome = genome,
            cpu = cpu,
            memory = memory
    }
    #step3
    call map_gap.Predata as Predata {
        input :
            workdir = workdir,
            sample = sample,
            hifi_fofn = hifi_fofn
    }
    call map_gap.map2gap as map2gap {
        input :
            hifi_fofn = Predata.new_hifidata,
            sample = sample,
            software_json = software_json,
            workdir = workdir,
            gap_genome = region.gap_genome,
            cpu = cpu,
            memory = memory
    }
    call map_gap.polish_gap as polish_gap {
        input :
            workdir = workdir,
            software_json = software_json,
            filt_bam = map2gap.filt_bam,
            filt_fa = map2gap.filt_fa,
            cpu = cpu,
            memory = memory
    }
    #step4
    call nextpolish1.nextPolish1 as nextPolish1 {
        input :
            software_json = software_json,
            workdir = polish_gap.wd2,
            sample = sample,
            lgs_fofn = lgs_fofn,
            sgs_fofn = sgs_fofn,
            hifi_fofn = hifi_fofn,
            genome = genome,
            cpu = cpu,
            memory = memory
    }
    #step5
    call nextpolish2fit.nextpolish2_fit as nextpolish2fit {
        input :
            genome = genome,
            workdir = nextPolish1.wd3,
            software_json = software_json,
            cpu = cpu,
            memory = memory
    }
    #step6
    call nextpolish2nofit.nextpolish2_nofit as nextpolish2nofit {
        input :
            software_json = software_json,
            workdir = nextpolish2fit.wd4,
            nomap_fa = region.nomap_fa,
            cpu = cpu,
            memory = memory
    }
    #step7
    call mergefasta.merge as merge {
        input :
            software_json = software_json,
            workdir = nextpolish2nofit.wd5,
            genome = genome,
            sgs_fofn = sgs_fofn,
            hifi_fofn = hifi_fofn,
            lgs_fofn = lgs_fofn,
            cpu = cpu,
            memory = memory
    }
    #step8
    call finallypolish.finally_polish as finallypolish {
        input :
            software_json = software_json,
            hifi_fofn = hifi_fofn,
            mergefasta = merge.mergefasta,
            workdir = workdir,
            cpu = cpu,
            memory = memory
    }
}
