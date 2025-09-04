version 1.0

# To run it: 
# module load caper
# caper run ExtractChromosome.wdl \
# -i /gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/ExtractChromosomeInput.json \
# -o /gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/ExtractChromosomeOptions.json \
# --backend-file /gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/backend.conf


workflow ExtractChromosomeWorkflow {
    input {
        File input_bam
        File input_bam_index
        String prefix
        String directory
    }

    call ExtractChromosome {
        input: input_bam = input_bam,
        input_bam_index = input_bam_index,
        prefix = prefix,
        directory = directory
    }
}

task ExtractChromosome {
    input {
        File input_bam
        File input_bam_index
        String prefix
        String directory
    }

    command {
        module load samtools/1.18
        samtools view -b -X ${input_bam} ${input_bam_index} 22 > ${directory}${prefix}_chr22_wdlplay.bam
        # samtools fixmate ${directory}${prefix}_chr22_wdlplay.bam ${directory}${prefix}_chr22_fixed_wdlplay.bam
        samtools view -f 2 -b ${directory}${prefix}_chr22_wdlplay.bam > ${directory}${prefix}_chr22_fixed_wdlplay.bam
        samtools sort ${directory}${prefix}_chr22_fixed_wdlplay.bam -o ${directory}${prefix}_chr22_fixed_sorted_wdlplay.bam
        samtools index ${directory}${prefix}_chr22_fixed_sorted_wdlplay.bam
    }

    runtime {
        memory: "10G"
        disks: "local-disk " + 200 + " HDD"
        cpu: 2
        preemptible: 1
    }

    output {
        File extracted_chromosome = "${directory}${prefix}_chr22_fixed_wdlplay.bam"
        File sorted_chromosome = "${directory}${prefix}_chr22_fixed_sorted_wdlplay.bam"
    }
}

# task IndexBam {
#     input {
#         File input_bam
#     }
    
#     command {
#         module load samtools
#         samtools index ${input_bam}
#     }
# }

# task SortBamByQname {
#     input {
#         File input_bam
#         String prefix
#     }

#     command {
#         module load samtools
#         samtools sort -n -b ${input_bam} ${prefix}_chr22_sorted.bam
#     }

#     output {
#         sorted_bam = "${prefix}_chr22_sorted.bam"
#     }
# }

# task ExtractChromosome {
#     input {
#         File input_bam
#         String prefix
#     }

#     command {
#         module load samtools
#         samtools view -f 2 -b ${input_bam} 22 > ${prefix}_chr22.bam
#     }

#     output {
#         extracted_bam = "${prefix}_chr22.bam"
#     }
# }
