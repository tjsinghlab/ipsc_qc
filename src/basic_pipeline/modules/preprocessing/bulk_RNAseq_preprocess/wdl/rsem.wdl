version 1.0

task rsem {

    input {
        File transcriptome_bam
        String rsem_reference
        String prefix

        String outdir

        Int memory=16
        Int disk_space
        Int num_threads=8
        Int num_preempt=0

        Int? max_frag_len
        String? estimate_rspd
        String? is_stranded
        String? paired_end
    }

    command {
        set -euo pipefail
        # mkdir rsem_reference
        # tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        mkdir -p ${outdir + "RSEM_outputs"}
        /src/run_RSEM.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            ${"--paired_end " + paired_end} \
            --threads ${num_threads} \
            -o ${outdir + "RSEM_outputs"} \
            ${rsem_reference} ${transcriptome_bam} ${prefix}
        gzip ${outdir + "RSEM_outputs"}/*.results
    }

    output {
        File genes="${outdir}RSEM_outputs/${prefix}.rsem.genes.results.gz"
        File isoforms="${outdir}RSEM_outputs/${prefix}.rsem.isoforms.results.gz"
    }

    runtime {
        singularity: "/gpfs/commons/groups/singh_lab/resources/RNAseq/gtex_rnaseq_V10.sif"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}