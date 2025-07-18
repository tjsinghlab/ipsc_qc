#!/bin/bash
#MYCOPLASMA DETECTION
#NEEDS FASTQ FILES AND MYCOPLASMA REFERENCE GENOME
#Use bowtie2 to align reads to mycoplasma genome, and assess alignment
#The alignment will be a proxy for level of contamination

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz
#--2025-01-24 16:08:05--  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz
#/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/GCF_000027325.1_ASM2732v1_genomic.fna.gz

module load bowtie2

#Build mycoplasma reference
bowtie2-build mycoplasma_genome.fna mycoplasma_index

#Align reads to mycoplasma reference
#Single-end:
bowtie2 -x mycoplasma_index -U sample_reads.fastq -S output.sam
#Paired-end:
bowtie2 -x mycoplasma_index -1 sample_R1.fastq -2 sample_R2.fastq -S output.sam

#Get alignment summary
bowtie2 -x mycoplasma_index -U sample_reads.fastq -S output.sam --no-unal 2> bowtie2_log.txt
cat bowtie2_log.txt
#Check alignment rate in the log to assess mycoplasma contamination

#Convert SAM to BAM to make our lives easier later
samtools view -bS output.sam > output.bam
samtools sort output.bam -o output_sorted.bam
samtools index output_sorted.bam