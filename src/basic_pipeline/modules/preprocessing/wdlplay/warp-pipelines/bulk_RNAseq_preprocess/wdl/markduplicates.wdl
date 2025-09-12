version 1.0

task markduplicates {

    input {
        File input_bam
        String prefix
        Int? max_records_in_ram
        Float? sorting_collection_size_ratio

        String outdir

        Float memory=16
        Int java_memory = floor(memory - 0.5)
        Int disk_space
        Int num_threads=8
        Int num_preempt=0

        String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")
    }

    command {
        set -euo pipefail
        python3 -u /src/run_MarkDuplicates.py ${input_bam} ${prefix} \
            -o ${outdir + "/Mark_duplicates_outputs"} \
            --memory ${java_memory} \
            ${"--max_records_in_ram " + max_records_in_ram} \
            ${"--sorting_collection_size_ratio " + sorting_collection_size_ratio}
        samtools index ${outdir + "/Mark_duplicates_outputs"}/${output_bam}
    }

    output {
        File bam_file = "${outdir}/Mark_duplicates_outputs/${output_bam}"
        File bam_index = "${outdir}/Mark_duplicates_outputs/${output_bam}.bai"
        File metrics = "${outdir}/Mark_duplicates_outputs/${prefix}.marked_dup_metrics.txt"
    }

    runtime {
        singularity: "/ref/gtex_rnaseq_V10.sif"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}
