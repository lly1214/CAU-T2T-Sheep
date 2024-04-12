version 1.0

import "tasks/map_gap.wdl" as map_gap

workflow Map_to_gap {
    input {
        File hifi_fofn
        String sample
        File software_json
        String workdir
        String gap_genome

        Int cpu = 20
        String memory = '50G'
    }
    call map_gap.predata as predata {
        input :
            workdir = workdir,
            sample = sample,
            hifi_fofn = hifi_fofn
    }
    call map_gap.map2gap as map2gap {
        input :
            hifi_fofn = predata.new_hifidata,
            sample = sample,
            software_json = software_json,
            workdir = workdir,
            gap_genome = gap_genome,
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
}