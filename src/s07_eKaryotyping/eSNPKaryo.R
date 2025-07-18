library(devtools)
library(eSNPKaryotyping)
library(zoo)
library(gplots)
library(patchwork)
library(ggplot2)

#BAM file/STAR outs location: "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/star_outs/"
#VCF file location: /gpfs/commons/groups/singh_lab/users/dongwang/RNA-variant-calling

# Before LOH analysis, the dbSNP files need to be edited using the Edit_dbSNP_Files function (done only once):
Edit_dbSNP_Files(Directory = "/Users/kjakubiak/Documents/chr/", File_Name = "chr", Organism = "Human")
# Argument: 1. Directory - the directory were the files are, one GTF file per chromosome
#           2. File_Name - the files name, without the number of the chromosomes
#           3. Organism - "Human" or "Mouse"

#EditVCF
Directory="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/"
Organism="Human"
print("Editing VCF File")
  Directory  
  Dir="/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/"
  file = "24246R-06-11.variant_filtered.vcf.gz"
  path = paste(Dir, file, sep="")
  readData = read.delim(path,as.is=T)
  readData=as.character(readData[-c(1:which(readData=="#CHROM")-1),1])
  
  jump = 10
  startChr = 1+jump
  startPos = 2+jump
  startInfo = 10+jump
  
  len = length(readData)
  chrRegex ="^chr(\\w+)$"
  infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"
  
  chrVector =  readData[startChr]
  posVector =  readData[startPos]
  infoVector = readData[startInfo]
  
  while (startInfo + jump < len) {
    startChr = startChr + jump
    startPos = startPos + jump
    startInfo = startInfo + jump
    chrVector = append(chrVector, readData[startChr])
    posVector = append(posVector, readData[startPos])
    infoVector = append(infoVector, readData[startInfo])
  }
  
  chrNum = gsub(chrRegex, "\\1", chrVector)
  if (Organism=="Human"){
    chrNum[chrNum=="X"]="23"
    chrNum[chrNum=="Y"]="24"}
  
  if (Organism=="Mouse"){
    chrNum[chrNum=="X"]="20"
    chrNum[chrNum=="Y"]="21"}
  
  chrNum =  as.numeric(chrNum)
  Karyotape = 10*abs(as.numeric(gsub(infoRegex, "\\1", infoVector))-as.numeric(gsub(infoRegex, "\\2", infoVector)))
  AD1 = as.numeric(gsub(infoRegex, "\\3", infoVector))
  AD2 = as.numeric(gsub(infoRegex, "\\4", infoVector))
  DP = as.numeric(gsub(infoRegex, "\\5", infoVector))
  
  posVector = as.numeric(posVector)
  
  table = data.frame("chr" = chrNum, "position" = posVector, "AD1" = AD1, "AD2" = AD2, "DP" = DP, "Karyotape" = Karyotape)
  
  fileName = "variantTable.csv"
  pathToSave = paste(Dir, fileName, sep="")
  table[is.na(table)] = 0
  write.table(table,pathToSave , sep="\t",row.names=F,quote=F)


#table<-read.csv("variantTable.csv")

# 5. Run the following lines:
table$chr=as.numeric(table$chr)
table=table[order(table$chr,table$position),]
table=table[table$chr>0,]

# 6. Run the MajorMinorCalc function
table2=MajorMinorCalc(Table = table,minDP = 20,maxDP = 1000000,minAF = 0.2)

PlotGenome(table2, Window=151, Organism="Human", Ylim=3, PValue=T)

#tbl=DeletionTable(Directory = "/Users/kjakubiak/Documents/variantTable.csv",Table = table2,dbSNP_Data_Directory = ,dbSNP_File_Name = ,Genome_Fa_dict = ,Organism = )

DeletionTable<-function(Directory,Table,dbSNP_Data_Directory,dbSNP_File_Name,Genome_Fa_dict,Organism){
  print("Reading SNPs table")
  i=1
  
  if(Organism=="Human"){mx=25}
  if(Organism=="Mouse"){mx=22}
  
  while (i <mx){
    print(paste("Chromosome:",i))
    chrTable = read.delim(paste(dbSNP_Data_Directory,dbSNP_File_Name,i,sep=""))
    if (i==1) {snpTable = chrTable}
    if (i>1) {snpTable = rbind(snpTable,chrTable)}
    i=i+1
  }
  
  
  table2=Table
  colnames(table2)[2]="start"
  print("Merging Tables")
  x=merge(snpTable,table2,by = c("chr","start"),all.x=T)
  x=x[order(x$chr,x$start),]
  
  
  setwd(Directory)
  
  print("Indexing BAM File")
  system("samtools index accepted_hits.bam")
  
  dict=read.csv(Genome_Fa_dict,as.is=T)
  dict_type=grep(pattern = "chr",x = dict[1,1])
  if(length(dict_type)==0){dict_type=0}
  
  
  i=1
  while (i<mx){
    print(paste("Chromosome ",i, "| ",Sys.time(),sep=""))
    x1=x[x$chr==i,]
    
    
    
    
    if(dict_type==0){
      loc=paste("chr",x1$chr[1],":",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")
      if(Organism=="Human"){
        if (i==23) {loc=paste("chrX:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==24) {loc=paste("chrY:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
      if(Organism=="Mouse"){
        if (i==20) {loc=paste("chrX:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==21) {loc=paste("chrY:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
    }
    
    if(dict_type==1){
      loc=paste(x1$chr[1],":",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")
      #loc=paste("chr", x1$chr[1], ":", x1$start[1], "-", x1$start[dim(x1)[1]], sep="")
      
      if(Organism=="Human"){
        if (i==23) {loc=paste("X:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==24) {loc=paste("Y:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
      if(Organism=="Mouse"){
        if (i==20) {loc=paste("X:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}
        if (i==21) {loc=paste("Y:",x1$start[1],"-",x1$start[dim(x1)[1]],sep="")}}
    }
    
    command=paste("samtools depth -r ",loc ," accepted_hits.bam > reads-per-position.txt",sep="")
    system(command)
    system ("awk -F \" \" '($3 >20){print $0}' reads-per-position.txt >reads-per-position2.txt")
    chr=read.delim("reads-per-position2.txt",header = F)
    colnames(chr)=c("chr","start","Depth")
    x2=merge(x1,chr,by="start")
    x2=cbind(x2,Depth_group=x2$Depth)
    x2$Depth_group[x2$Depth<50] ="20-50"
    x2$Depth_group[x2$Depth>49 & x2$Depth<100] ="50-100"
    x2$Depth_group[x2$Depth>99 & x2$Depth<200] ="100-200"
    x2$Depth_group[x2$Depth>199 & x2$Depth<500] ="200-500"
    x2$Depth_group[x2$Depth>499] =">500"
    x2$Depth_group=as.factor(x2$Depth_group)
    if (i==1) {tbl=x2}
    if (i>1) {tbl=rbind(tbl,x2)}
    i=i+1
  }
  
  
  print("Writing Table")
  tbl[is.na(tbl)]=0
  write.table(tbl,"Deletions.txt" , sep="\t",row.names=F,quote=F)
  return(tbl)
}

# 8. intersect the observed SNPs with the common SNPs table from dbSNPs, Creates file with the LOH data called Deletions.txt
#this requires the accepted_hits.bam output from tophat, so I've renamed the aligned.out.bam from STAR to accepted_hits.bam to match what the function is expecting
tbl=DeletionTable(Directory = "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/tmp/",Table = table2,dbSNP_Data_Directory = "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/chr/",dbSNP_File_Name = "Edited_Common_chr",Genome_Fa_dict = "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ekaryo/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna",Organism = "Human")
# Argument: 1. Directory - The Path of to the BAM Directory, containing the variantTable.csv file
#           2. Table - The variable containing the output of the MajorMinorCalc function
#           3. dbSNP_Data_Directory - The path for the directory where the edited dbSNP file are (created previously by the Edit_dbSNP_Files function)
#           4. dbSNP_File_Name - The edited dbSNP file names, without the chromosome number
#           5. Genome_Fa_dict - the path for the dictionary file created for the whole genome FASTA file.
#           6. Organism - "Human" or "Mouse"


# To run the LOH analysis again, just read the Deletions.txt file into a variable using the read.delim function


Plot_Zygosity_Sinle(Table = tbl,Organism = "Human")

Plot_Zygosity_Blocks(Table = tbl,Window = 1500000,Max = 6,Max2 = 60,Organism = "Human")

