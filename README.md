# iPSC QC Pipeline

<img width="2327" height="380" alt="Screenshot 2025-10-09 at 2 42 57 PM" src="https://github.com/user-attachments/assets/61f7b78d-afed-43fe-9e13-06de2e3377c6" />

This pipeline includes: 
- Processing raw fastq files (bulk RNA sequencing) using fastqc, STAR alignment, and RSEM
- GATK germline variant calling
- Cancer mutation calling (in select oncogenes)
- eSNPKaryotyping
- Mycoplasma detection (*Mycoplasma fermentans* and *Mycoplasma orale*)
- PACNet classification
- Outlier assessment
  
This pipeline expects a directory of fastq.gz (paired end: R1 and R2) files from a bulk RNA sequencing run.

See [wiki](https://github.com/tjsinghlab/ipsc_qc/wiki/Running-the-Pipeline) for more details.

---

## Table of Contents

- [Running the Pipeline](#running-the-pipeline)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Directory Structure](#directory-structure)
- [Notes](#notes)
- [Sources and Citations](#citations)

---

## Running the Pipeline
### Example Command

```
#Download docker image as singularity image (for HPC compatibility):

singularity pull ipsc_qc.sif docker://kjakubiak16/ipsc_qc:latest

#Run the singularity command. First bind the directories, then supply bound directories as arguments.

singularity exec \
  -B ./path/to/fastqs:/data \
  -B ./path/to/ref_dir:/ref \
  -B ./path/to/COSMIC_DB:/cosmic \
  -B ./path/to/output_dir:/output \
  ipsc_qc.sif \
  /pipeline/pipeline_runner.sh \
    --fastq_dir /data \
    --output_dir /output \
    --ref_dir /ref \
    --cosmic_dir /cosmic
```

## Requirements
### Reference COSMIC Database Files (provide in --cosmic_dir argument; default is ./cosmic)
#### Must be provided by user
- `Cosmic_CancerGeneCensus_Tsv_v101_GRCh38.tar`  
  *Cosmic Database; MUST BE DOWNLOADED BY USER AND PROVIDED IN COSMIC_DIR ARGUMENT. If not, cancer mutation calling will not be performed.*

### Reference Files (place in `ref/` directory; default is ./ref)
- `chr/`
  *Directory with dbSNP build 142 common SNP files (GTF format), one per chromosome. Human: X = 23, Y = 24. Source: UCSC Table Browser, snp142common table.*
- `genes.txt`
  *File containing a list of oncogenes for assessment. Lives in this repo; download into your /ref dir (the one you provide to the --ref_dir argument), and add whichever genes you want analyzed against COSMIC.*
- GTEX and RSEM References (`Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta`, `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai`, `Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict`, `gencode.v39.GRCh38.genes.collapsed_only`  
  *(requester pays enabled; cannot be universally downloaded within pipeline))*

#### If not included in user-provided reference directory (ref_dir), the following files will be automatically downloaded to ref_dir:
- `star_index_oh75`
  *reference files for STAR*
- `rsem_reference`
  *refDict*
- `Homo_sapiens_assembly38.dbsnp138.vcf`
  *dbSnpVcf*
- `Homo_sapiens_assembly38.dbsnp138.vcf.idx` 
  *dbSnpVcfIndex*
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
  *knownVcf 1*
- `Homo_sapiens_assembly38.known_indels.vcf.gz`
  *knownVcf 2*
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi`
  *knownVcfIndex 1*
- `Homo_sapiens_assembly38.known_indels.vcf.gz.tbi`
  *knownVcfIndex 2*
- `gencode.v39.GRCh38.genes.collapsed_only.gtf`
  *annotationsGTF*
- `images/gatk_4.6.1.0.sif`
  *gatk4_docker image*
- `GCF_000420105.1_ASM42010v1_genomic.fna.gz`  
  *Mycoplasma orale genome; downloaded from `wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/420/105/GCF_000420105.1_ASM42010v1/GCF_000420105.1_ASM42010v1_genomic.fna.gz`*
- `GCF_003704055.1_ASM370405v1_genomic.fna.gz`
  *Mycoplasma fermentans genome; downloaded from `wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/055/GCF_003704055.1_ASM370405v1/GCF_003704055.1_ASM370405v1_genomic.fna.gz`*
- `Hs_expTrain_Jun-20-2017.rda`  
  *PACNet training data; downloaded from `aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda . --no-sign-request`*
- `Hs_stTrain_Jun-20-2017.rda`  
  *PACNet training data; downloaded from `aws s3 cp s3://cellnet-rnaseq/ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda . --no-sign-request`*
- `GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna`  
  *NCBI reference genome FASTA file; downloaded from `ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38.p13`*

### Tools
The pipeline depends on the following command line tools:
#### Must be installed by user (Ensure these tools are installed or available via a module system (ex. HPC users)):
- `singularity`
- `parallel`
- `squashfuse`

## Inputs

| Argument         | Description                                                                | Default                        |
|------------------|----------------------------------------------------------------------------|--------------------------------|
| `--ref_dir`      | Directory for reference files                                              | `/refs`                        |
| `--fastq_dir`    | Directory containing input files (fastq.gz sequencing files) for a run     | `./fastq`                      |
| `--output_dir`   | Path to desired output directory                                           | `./outputs`                    |
| `--cosmic_dir`   | Directory containing downloaded COSMIC database files.                     | `cosmic`                       |
| `--project`      | Project name                                                               | `default_project`              |

**Note:**  
File base names should be unique for each sample and consistent across file types  
(e.g., `sample1.vcf`, `sample1.bam`, `sample1.genes.results`, `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`)

---

## Outputs
### Output Directory Structure

```
output_dir
  ├── logs
  │    ├── ekaryo_sample_01.log
  │    ├── ekaryo_sample_02.log
  │    ├── cancer_mutation_mapping_sample_01.log
  │    ├── cancer_mutation_mapping_sample_02.log
  │    ├── html_summary_sample_01.log
  │    ├── html_summary_sample_02.log
  │    ├── pacnet_all.log
  │    ├── outliers_sample_01.log
  │    ├── outliers_sample_02.log
  │    ├── sample_01_wdl2.log
  │    ├── sample_02_wdl2.log
  │    ├── sample_01_bulk_preprocess.log
  │    └── sample_02_bulk_preprocess.log  
  ├── pacnet
  │    ├── PACNet_heatmap.png
  │    ├── classifier.rda
  │    ├── query_meta.csv
  │    ├── query_matrix.csv
  │    ├── metrics.json
  │    ├── classification_validation_hm.png
  │    └── classification_scores.csv
  ├── outlier_analysis
  │    └── PCA_pacnet_scores.pdf
  ├── sample_01
  │    ├── mycoplasma
  │    │   ├── sample_01_flagstat.txt
  │    │   ├── sample_01_sorted.bam.bai
  │    │   ├── sample_01_sorted.bam
  │    │   ├── sample_01_bowtie2.log
  │    │   ├── mycoplasma_alignment_summary.png
  │    │   ├── mycoplasma_alignment_summary.pdf
  │    │   └── mycoplasma_alignment_stats.tsv
  │    ├── eSNPKaryotyping
  │    │   ├── sample_01_PlotGenome.png
  │    │   └── sample_01_variantTable.csv
  │    ├── cosmic_calling
  │    │   ├── CancerMutationPlot.pdf
  │    │   └── sample_01_CancerMutations.tsv
  │    ├── variant_calling
  │    │   ├── star_out
  │    │   │   ├── sample_01.Aligned.toTranscriptome.out.bam
  │    │   │   ├── sample_01.Aligned.out.bam
  │    │   │   ├── sample_01.Aligned.sortedByCoord.out.bam
  │    │   │   ├── sample_01.Aligned.sortedByCoord.out.bam.bai
  │    │   │   ├── sample_01.ReadsPerGene.out.tab.gz
  │    │   │   ├── sample_01.Chimeric.out.sorted.bam
  │    │   │   └── sample_01.Chimeric.out.sorted.bam.bai
  │    │   ├── RSEM_outputs
  │    │   │   ├── sample_01.rsem.isoforms.results.gz
  │    │   │   ├── sample_01.rsem.genes.results.gz
  │    │   │   └── sample_01.rsem.stat  
  │    │   ├── QC_outputs
  │    │   │   ├── sample_01_R1.fastqc_data.txt.gz
  │    │   │   ├── sample_01_R2.fastqc_data.txt.gz
  │    │   │   ├── sample_01_R1_fastqc.zip
  │    │   │   ├── sample_01_R2_fastqc.zip
  │    │   │   ├── sample_01_R1_fastqc.html
  │    │   │   └── sample_01_R2_fastqc.html
  │    │   ├── fastqc_out
  │    │   │   ├── sample_01_R1.fastqc_data.txt.gz
  │    │   │   ├── sample_01_R2.fastqc_data.txt.gz
  │    │   │   ├── sample_01_R1_fastqc.zip
  │    │   │   ├── sample_01_R2_fastqc.zip
  │    │   │   ├── sample_01_R1_fastqc.html
  │    │   │   └── sample_01_R2_fastqc.html
  │    │   └── Mark_duplicates_outputs
  │    │   │   ├── sample_01.Aligned.sortedByCoord.out.md.bam.bai
  │    │   │   ├── sample_01.marked_dup_metrics.txt
  │    │   │   └── sample_01.Aligned.sortedByCoord.out.md.bam
  │    ├── preprocessing
  │    ├── logs
  │    ├── variant_calling
  │    │   ├── sample_01.variant_filtered.vcf.gz
  │    │   ├── sample_01.variant_filtered.vcf.gz.tbi
  │    │   ├── sample_01.aligned.duplicates_marked.recalibrated.bam
  │    │   └── sample_01.aligned.duplicates_marked.recalibrated.bam.bai
  └── sample_02
  


```

## Directory structure
### Basic pipeline repo structure
```
src
  ├── basic_pipeline
       ├── docker_image_structure
       ├── docker-compose.yaml
       ├── Dockerfile_slim
       ├── modules
       │   ├── cancer_mutation_calling
       │   │   └── COSMIC_cancer_mutation_calling.r
       │   ├── eSNPKaryotyping
       │   │   └── run_eSNPKaryotyping.R
       │   ├── mycoplasma_detection
       │   │   └── detect_mycoplasma.sh
       │   ├── outlier_detection
       │   │   └── outlier_detection.R
       │   ├── PACNet
       │   │   └── run_pacnet.R
       │   └── preprocessing
       │       └── wdlplay
       │           ├── LICENSE.md
       │           ├── Makefile
       │           ├── MANIFEST.in
       │           ├── README.md
       │           ├── requirements.txt
       │           ├── settings.json
       │           ├── setup.cfg
       │           ├── setup.py
       │           ├── versioneer.py
       │           ├── warp-pipelines
       │           │   ├── bulk_RNAseq_preprocess
       │           │   │   ├── RNAseq_pipeline_fastq.wdl
       │           │   │   ├── run_wdl.py
       │           │   │   └── wdl
       │           │   │       ├── fastqc.wdl
       │           │   │       ├── markduplicates.wdl
       │           │   │       ├── rnaseqc2.wdl
       │           │   │       ├── rsem.wdl
       │           │   │       └── star.wdl
       │           │   └── GATK_variant_calling
       │           │       ├── gatk4-rna-best-practices.wdl
       │           │       ├── gatk4-rna-germline-calling_run.py
       │           │       └── gatk4-rna-germline-variant-calling.inputs.json
       │           └── wdlplay
       │               ├── __init__.py
       │               ├── _version.py
       │               ├── db.json
       │               └── wdlplayer.py
       ├── pipeline_runner.sh
       ├── README.md
       ├── ref_files
       │   └── genes.txt
       ├── report_builder.R
       └── tarballs
           └── CellNet_master.tar.gz
```
---

## Notes
- This pipeline is memory-intensive. Implementation on HPC is recommended.
- WDL scripts will utilize cromwell, but a caper server will not be started. Backend is local. All packages and scripts are contained within the image; user does not need to install anything beyond what is included in the [requirements](#requirements) section.

---

## Citations

- **COSMIC: the Catalogue Of Somatic Mutations In Cancer**: Tate JG, Bamford S, Jubb HC, Sondka Z, Beare DM, Bindal N, Boutselakis H, Cole CG, Creatore C, Dawson E, Fish P, Harsha B, Hathaway C, Jupe SC, Kok CY, Noble K, Ponting L, Ramshaw CC, Rye CE, Speedy HE, Stefancsik R, Thompson SL, Wang S, Ward S, Campbell PJ, Forbes SA. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res. 2019 Jan 8;47(D1):D941-D947. doi: 10.1093/nar/gky1015. PMID: 30371878; PMCID: PMC6323903.
- **PACNet (Platform-Agnostic CellNet)**: Lo EKW, Velazquez JJ, Peng D, Kwon C, Ebrahimkhani MR, Cahan P. Platform-agnostic CellNet enables cross-study analysis of cell fate engineering protocols. Stem Cell Reports. 2023 Aug 8;18(8):1721-1742. doi: 10.1016/j.stemcr.2023.06.008. Epub 2023 Jul 20. PMID: 37478860; PMCID: PMC10444577.
- **eSNPKaryotyping**: Weissbein U, Schachter M, Egli D, Benvenisty N. Analysis of chromosomal aberrations and recombination by allelic bias in RNA-Seq. Nat Commun. 2016 Jul 7;7:12144. doi: 10.1038/ncomms12144. PMID: 27385103; PMCID: PMC4941052.
- **GATK, Docker, WDL**: Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
- **fastqc**: Li, B., Dewey, C.N. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12, 323 (2011). https://doi.org/10.1186/1471-2105-12-323
- **RSEM**: Li, B., Dewey, C.N. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12, 323 (2011). https://doi.org/10.1186/1471-2105-12-323
- **STAR Alignment**: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.
- **parallel**: O. Tange (2018): GNU Parallel 2018, March 2018, https://doi.org/10.5281/zenodo.1146014.
- **Bowtie**: Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923
