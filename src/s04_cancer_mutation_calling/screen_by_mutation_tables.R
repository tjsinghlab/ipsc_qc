library(data.table)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)
setwd("/gpfs/commons/groups/singh_lab/users/kjakubiak")
VCF_FILE="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/24246R-06-15.variant_filtered.vcf.gz"
VCF_FILE="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/24246R-06-17.variant_filtered.vcf.gz"
cos="/gpfs/commons/groups/singh_lab/users/kjakubiak/COSMIC_DB/Cosmic_MutantCensus_v101_GRCh38.tsv.gz"
GENE_LIST="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/"

tmp_vcf<-readLines(VCF_FILE)
tmp_vcf_data<-read.table(VCF_FILE, stringsAsFactors = FALSE)

# filter for the columns names
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
names(tmp_vcf_data)<-vcf_names

library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

brca2<-read_excel_allsheets("ipsc/mutation_screen/published_tables/BRCA2_Huang_Couch.xlsx")
#brca2_3<-brca2$`Table S3`
#brca2_6b<-brca2$`Table S6B`
brca2<-brca2$`Table S6B`

#colnames(brca2_3)<-brca2_3[1,]
#colnames(brca2_6b)<-brca2_6b[1,]
colnames(brca2)<-brca2[1,]

#brca2_3<-brca2_3[-c(1),]
#brca2_6b<-brca2_6b[-c(1),]
brca2<-brca2[-c(1),]

#brca2_3<-brca2_3[-c(1,3:4,11:21)]
#brca2_6b<-brca2_6b[-c(1,3:4,11:21)]
brca2<-brca2[-c(1,3,11:21)]
colnames(brca2)[1]<-"CHROMOSOME"
brca2$POS<-str_split_i(brca2$`Annotation(GRCh38)`, "_", 2)
#colnames(brca2)[3]<-"POS"

brca2$REF<-str_split_i(brca2$`Annotation(GRCh38)`, "_", 3)
brca2$ALT<-str_split_i(brca2$`Annotation(GRCh38)`, "_", 4)

tp53<-read_excel_allsheets("ipsc/mutation_screen/published_tables/TP53_Merkle_Eggan.xlsx")
tp53<-tp53$'1_P53_mutations_summary'

tp53$CHROMOSOME<-str_split_i(tp53$`Genomic position (hg19)`, ":", 1)
tp53$POS<-str_split_i(tp53$`Genomic position (hg19)`, ":", 2)
tp53<-tp53[-c(1,2,4:8,12:24)]

tp53$REF<-str_split_i(tp53$`Genomic DNA alteration`, ">", 1)
tp53$ALT<-str_split_i(tp53$`Genomic DNA alteration`, ">", 2)

cosmic = fread("/gpfs/commons/groups/singh_lab/users/kjakubiak/COSMIC_DB/Cosmic_MutantCensus_v101_GRCh38.tsv.gz")
cosmic$POS<-cosmic$GENOME_START
cosmic$REF<-cosmic$GENOMIC_WT_ALLELE
cosmic$ALT<-cosmic$GENOMIC_MUT_ALLELE

VCF_FILE="/gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling/24246R-06-11.variant_filtered.vcf.gz"
tmp_vcf<-readLines(VCF_FILE)
tmp_vcf_data<-read.table(VCF_FILE, stringsAsFactors = FALSE)

# filter for the columns names
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
names(tmp_vcf_data)<-vcf_names
vcf<-tmp_vcf_data
colnames(vcf)[1]<-"CHROMOSOME"
vcf$CHROMOSOME<-str_split_i(vcf$CHROMOSOME, "r", 2)
vcf$CHROMOSOME<-as.numeric(vcf$CHROMOSOME)

merge_tp53_2<-merge(vcf, tp53, by = c("CHROMOSOME", "POS", "REF", "ALT"))
merge_brca2_2<-merge(vcf, brca2, by = c("CHROMOSOME", "POS", "REF", "ALT"))
merge_cosmic_2<-merge(vcf, cosmic, by = c("CHROMOSOME", "POS", "REF", "ALT"))

