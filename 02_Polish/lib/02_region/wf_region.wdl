version 1.0

import "tasks/region.wdl" as region

workflow Region {
    input {
        File hifi_fofn  
        String sample 
        File software_json 
        String workdir 
        File genome 

        Int cpu = 20
        String memory = '50G'
    }
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
            workdir = workdir,
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
}
