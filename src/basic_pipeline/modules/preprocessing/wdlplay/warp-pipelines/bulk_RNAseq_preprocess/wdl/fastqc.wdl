version 1.0

task fastqc {

    input {
        File fastq1

        String outdir

        Float memory=16
        Int disk_space
        Int num_threads=8
        Int num_preempt=0

        String fastq1_name = sub(sub(basename(fastq1), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )
    }

    command {
        set -euo pipefail
        mkdir -p ${outdir + "/fastqc_out"}
        fastqc ${fastq1} \
            --threads ${num_threads} \
            --outdir ${outdir + "/fastqc_out"}
        unzip -p ${outdir + "/fastqc_out"}/${fastq1_name}_fastqc.zip ${fastq1_name}_fastqc/fastqc_data.txt | gzip > ${outdir + "/fastqc_out"}/${fastq1_name}.fastqc_data.txt.gz
    }

    output {
        File fastq1_fastqc_html = "${outdir}/fastqc_out/${fastq1_name}_fastqc.html"
        File fastq1_fastqc_zip =  "${outdir}/fastqc_out/${fastq1_name}_fastqc.zip"
        File fastq1_fastqc_data = "${outdir}/fastqc_out/${fastq1_name}.fastqc_data.txt.gz"
    }

    runtime {
        singularity: "/pipeline/modules/gtex_rnaseq_V10.sif"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}