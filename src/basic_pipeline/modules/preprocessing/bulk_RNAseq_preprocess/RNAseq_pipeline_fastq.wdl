version 1.0

import "wdl/fastqc.wdl" as fastqc_wdl
import "wdl/star.wdl" as star_wdl
import "wdl/markduplicates.wdl" as markduplicates_wdl
import "wdl/rsem.wdl" as rsem_wdl
import "wdl/rnaseqc2.wdl" as rnaseqc_wdl

workflow rnaseq_pipeline_fastq_workflow {
    input {
        Array[File] fastqs
        File fastq1=fastqs[0]
        File fastq2=fastqs[1]
        String prefix=basename(fastq1, "_R1.fastq.gz")
        String star_index_oh75
        String rsem_reference
        File genes_gtf
        String outdir
    }

    Int disk_space=ceil(10 * size(fastq1, "GiB"))
    Int num_threads = 8
    Int num_preempt = 0

    call fastqc_wdl.fastqc as fastqc {
        input: fastq1=fastq1,
        fastq2=fastq2,
        outdir=outdir,
        disk_space=disk_space
    }

    call star_wdl.star as star {
        input: fastq1=fastq1, 
        fastq2=fastq2, 
        prefix=prefix,
        star_index=star_index_oh75,
        outdir=outdir
    }

    call markduplicates_wdl.markduplicates as markduplicates {
        input: input_bam=star.bam_file, 
        prefix=prefix,
        disk_space=disk_space,
        outdir=outdir
    }

    call rsem_wdl.rsem as rsem {
        input: transcriptome_bam=star.transcriptome_bam, 
        prefix=prefix, 
        rsem_reference=rsem_reference,
        disk_space=disk_space,
        outdir=outdir
    }

    call rnaseqc_wdl.rnaseqc2 as qc {
        input: bam_file=markduplicates.bam_file, 
        sample_id=prefix, 
        genes_gtf=genes_gtf,
        disk_space=disk_space,
        outdir=outdir
    }

    output {
        File fastq1_fastqc_html = fastqc.fastq1_fastqc_html
        File fastq1_fastqc_zip =  fastqc.fastq1_fastqc_zip
        File fastq1_fastqc_data = fastqc.fastq1_fastqc_data
        File fastq2_fastqc_html = fastqc.fastq2_fastqc_html
        File fastq2_fastqc_zip =  fastqc.fastq2_fastqc_zip
        File fastq2_fastqc_data = fastqc.fastq2_fastqc_data
        File star_bam_file = star.bam_file
        File star_bam_index = star.bam_index
        File star_transcriptome_bam = star.transcriptome_bam
        File star_chimeric_junctions = star.chimeric_junctions
        File star_chimeric_bam_file = star.chimeric_bam_file
        File star_chimeric_bam_index = star.chimeric_bam_index
        File star_read_counts = star.read_counts
        File star_junctions = star.junctions
        File star_junctions_pass1 = star.junctions_pass1
        Array[File] star_logs = star.logs

        File md_bam_file = markduplicates.bam_file
        File md_bam_index = markduplicates.bam_index
        File md_metrics = markduplicates.metrics

        File rsem_genes = rsem.genes
        File rsem_isoforms = rsem.isoforms

        File qc_gene_tpm = qc.gene_tpm
        File qc_gene_counts = qc.gene_counts
        File qc_exon_counts = qc.exon_counts
        File qc_metrics = qc.metrics
        File qc_gc_content = qc.gc_content
        File qc_insertsize_distr = qc.insertsize_distr
    }
}