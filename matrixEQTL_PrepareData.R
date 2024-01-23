@ -0,0 +1,161 @@
################################################
##   The following files should be in the same directory
#  paste0(FilePrefix,".fam")
#  paste0(FilePrefix,".eigenvec")
#  paste0(FilePrefix,".expression.rds")     expression matrix
#  paste0(FilePrefix,".expression.csv")     Phenotype data
#  paste0(FilePrefix,".chr",i,".raw")       Genotype data for each chromosome
#  paste0(FilePrefix,".GeneLocation.csv")   Location of genes (geneid,chr,start,end)
####################################################################################

args<-commandArgs(trailingOnly = TRUE)

suppressMessages(library(stringr))
suppressMessages(library(data.table))

InDir<-args[1] 
OutDir <- args[2] 
covar_fact <- args[3] 
covar_num <- args[4] 
FilePrefix <- args[5] 
chr <- args[6]


no.cov=F

if(chr=="all"){
  chr <- seq(1,22,1)
}else{
  chr <- as.numeric(str_split(args[6], pattern = ",", simplify = T)[1,])
}

samples <- read.table(file=paste0(InDir,"/",FilePrefix,".fam"),sep=" ",header=F,stringsAsFactors=F)
rownames(samples) <- samples$V2

print("Reading genotype eigenvectors...")
eigenvec <- read.table(file=paste0(OutDir,"/",FilePrefix,".eigenvec"),header=F,row.names = 2)[,c(-1)]

print("Reading expression data...")
exp_all <- as.data.frame(readRDS(paste0(InDir,"/",FilePrefix,".expression.rds")))

if(!all(sapply(list(colnames(exp_all),rownames(eigenvec)), FUN = identical, rownames(samples)))){
  shared_names <- Reduce(intersect, list(colnames(exp_all),rownames(samples),rownames(eigenvec)))
  exp_all <- exp_all[,colnames(exp_all) %in% shared_names]
  samples <- samples[rownames(samples) %in% shared_names,]
  
  exp_all <- exp_all[shared_names]
  samples <- samples[match(shared_names,rownames(samples)),]
  eigenvec <- eigenvec[match(shared_names,rownames(eigenvec)),]
}

if(file.exists(paste0(InDir,"/",FilePrefix,".expression.csv"))){
  exp_sample <- read.csv(paste0(InDir,"/",FilePrefix,".expression.csv"),row.names=1,stringsAsFactors = F) 
  
  if(!identical(rownames(exp_sample) , colnames(exp_all))){
    shared_names <- intersect(colnames(exp_all),rownames(exp_sample))
    exp_sample <- exp_sample[rownames(exp_sample) %in% shared_names,]
    exp_sample <- exp_sample[match(shared_names,rownames(exp_sample)),]
  }
}else{
  no.cov=T
}
snps_all=vector(mode = "list",length = 22)
snps_loc=vector(mode = "list",length = 22)
for (i in chr) {
  print(paste("Reading Genotype data chr",i,"..."))
  if(!file.exists(paste0(OutDir,"/",FilePrefix,".chr",i,".genotype.txt"))){
  
    snps_all[[i]]<-fread(paste0(OutDir,"/",FilePrefix,".chr",i,".raw"), data.table = FALSE,nThread=16)
    
    snps_all[[i]] <- t(snps_all[[i]])
    snps_all[[i]] <- cbind.data.frame(rownames(snps_all[[i]]),snps_all[[i]])
    names(snps_all[[i]]) <- c("snpid",snps_all[[i]][2,c(-1)])
    snps_all[[i]] <- snps_all[[i]][c(-1:-6),]
    write.table(snps_all[[i]],file=paste0(OutDir,"/",FilePrefix,".chr",i,".genotype.txt"),sep="\t",row.names = F,col.names = T,quote = F)
  } else {
    snps_all[[i]]<-fread(paste0(OutDir,"/",FilePrefix,".chr",i,".genotype.txt"), data.table = FALSE,nThread=16)
  }
  if(!identical(colnames(exp_all),colnames(snps_all[[i]])[-1])){
    index <- match(colnames(exp_all),colnames(snps_all[[i]])[-1])
    snps_all[[i]] <- snps_all[[i]][,c(1,index+1)]
  }
  snps_loc[[i]] <- str_split(snps_all[[i]]$snpid,pattern=':',simplify=T)[,c(1,2)]
  snps_loc[[i]] <- cbind.data.frame(snps_all[[i]]$snpid,snps_loc[[i]])
  names(snps_loc[[i]]) <- c("snpid","chr","pos")
} 

print("Reading Gene location data...")
gene_loc_all <- read.csv(file=paste0(InDir,"/GeneLocation.csv"),stringsAsFactors = F,row.names = 1)
#names(gene_loc_all) <- tolower(names(gene_loc_all))
gene_loc_all <- cbind.data.frame(rownames(gene_loc_all),gene_loc_all$chr,gene_loc_all$start,gene_loc_all$end)
names(gene_loc_all) <- c("geneid","chr","left","right")

if(!identical(row.names(exp_all) , gene_loc_all$geneid)){
  index <- match(row.names(exp_all) , gene_loc_all$geneid)
  gene_loc_all <- gene_loc_all[index,]
}

exp_all <- cbind.data.frame(gene_loc_all$geneid,exp_all)
names(exp_all)[1] <- "geneid"

print(paste("Are the gene ids matched in gene expression and location data?"))
print(ifelse(identical(exp_all$geneid,gene_loc_all$geneid),"Yes","NO"))

print("Saving Expression data...")
write.table(exp_all,file = paste0(OutDir,"/",FilePrefix,".exp.txt"),quote = F,col.names = T,row.names = F,sep = '\t')

print("Saving gene location data...")
write.table(gene_loc_all,file = paste0(OutDir,"/",FilePrefix,".gene.loc.txt"),quote = F,col.names = T,row.names = F,sep = '\t')

for (i in chr) {
  
  print(paste("Are the SNP ids matched in in chr",i," data?"))
  print(ifelse(identical(snps_all[[i]]$snpid,snps_loc[[i]]$snpid),"Yes","NO"))
  
  print(paste("Are the sample ids matched in expression data and chr",i," SNP data?"))
  print(ifelse(identical(colnames(exp_all)[-1], colnames(snps_all[[i]])[-1]),"Yes","NO"))
  
  print(paste("Saving Chr",i,"..."))
  write.table(snps_all[[i]],file = paste0(OutDir,"/",FilePrefix,".snps",".chr",i,".txt"),quote = F,col.names = T,row.names = F,sep = '\t')
  write.table(snps_loc[[i]],file = paste0(OutDir,"/",FilePrefix,".snps.loc",".chr",i,".txt"),quote = F,col.names = T,row.names = F,sep = '\t')

}

covar_fact = str_split(covar_fact,pattern = ',',simplify = T)[1,]
covar_num = str_split(covar_num,pattern = ',',simplify = T)[1,]

if((covar_fact!="")&(covar_num!="")){
  covariates <- cbind.data.frame(exp_sample[,covar_fact],exp_sample[,covar_num],eigenvec$V3,eigenvec$V4,eigenvec$V5)
  names(covariates) <- c(covar_fact,covar_num,"PC1","PC2","PC3")
}else{
  if((covar_fact!="")){
    covariates <- cbind.data.frame(exp_sample[,covar_fact],eigenvec$V3,eigenvec$V4,eigenvec$V5)
    names(covariates) <- c(covar_fact,"PC1","PC2","PC3")
  }else{
    if((covar_num!="")){
      covariates <- cbind.data.frame(exp_sample[,covar_num],eigenvec$V3,eigenvec$V4,eigenvec$V5)
      names(covariates) <- c(covar_num,"PC1","PC2","PC3")
    }else{
      no.cov=T
    }
  }
}

if(!no.cov){
  print("Generating covariates file...")
  rownames(covariates) <- rownames(exp_sample)
  for (c in covar_fact) {
    covariates[,c] = as.numeric(as.factor(covariates[,c]))
  }
  covariates <- t(covariates)
  covariates <- cbind.data.frame(rownames(covariates),covariates)
  names(covariates)[1] <- "id"
  
  print("Sample ids matched in covariate file with the others?")
  print(ifelse(all(sapply(list(colnames(exp_all)[-1], rownames(exp_sample),rownames(eigenvec),colnames(covariates)[-1]), FUN = identical, rownames(samples))),"Yes","NO"))
  #True
  
  print("Saving covaraites file...")
  write.table(covariates,file = paste0(OutDir,"/",FilePrefix,".covariates.txt"),quote = F,col.names = T,row.names = F,sep = '\t')
}
print("All done!")
