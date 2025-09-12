version 1.0

task rnaseqc2 {

    input {
        File bam_file
        File genes_gtf
        String sample_id
        String? strandedness
        File? intervals_bed
        File? reference_fasta
        File? reference_fasta_index
        String? flags

        String outdir

        Int memory=16
        Int disk_space
        Int num_threads=8
        Int num_preempt=0
    }

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
        mkdir -p ${outdir + "QC_outputs"}
        touch ${outdir + "QC_outputs"}/${sample_id}.fragmentSizes.txt
        touch ${outdir + "QC_outputs"}/${sample_id}.gc_content.tsv
        rnaseqc ${genes_gtf} ${bam_file} ${outdir + "QC_outputs"} -s ${sample_id} ${"--bed " + intervals_bed} ${"--stranded " + strandedness} ${"--fasta " + reference_fasta} -vv ${flags} 
        echo "  * compressing outputs"
        gzip ${outdir + "QC_outputs"}/*.gct
        echo $(date +"[%b %d %H:%M:%S] done")
    }

    output {
        File gene_tpm = "${outdir}/QC_outputs/${sample_id}.gene_tpm.gct.gz"
        File gene_counts = "${outdir}/QC_outputs/${sample_id}.gene_reads.gct.gz"
        File exon_counts = "${outdir}/QC_outputs/${sample_id}.exon_reads.gct.gz"
        File metrics = "${outdir}/QC_outputs/${sample_id}.metrics.tsv"
        File gc_content = "${outdir}/QC_outputs/${sample_id}.gc_content.tsv"
        File insertsize_distr = "${outdir}/QC_outputs/${sample_id}.fragmentSizes.txt"
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