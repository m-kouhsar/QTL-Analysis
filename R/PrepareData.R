args<-commandArgs(trailingOnly = TRUE)

suppressMessages(library(stringr))
suppressMessages(library(data.table))

fam.file <- trimws(args[1])
eigenvec.file <- trimws(args[2])
exp.rds.file <- trimws(args[3])
exp.pheno.file <- trimws(args[4])
geneLocation.csv.file <- trimws(args[5])
covar.fact <- trimws(args[6])
covar.num <- trimws(args[7])
OutPrefix <- trimws(args[8]) 
chr <- trimws(args[9])
overwrite <- ifelse(trimws(tolower(args[10]))=="yes" , T ,F)

cat("############################### PrepareData.R script ###########################################\n")
cat("Input arguments:\n")
cat("     fam file: ",fam.file,"\n")
cat("     eigenvec file: ",eigenvec.file,"\n")
cat("     Methylation/Expression matrix file: ",exp.rds.file,"\n")
cat("     Phenotype file: ",exp.pheno.file,"\n")
cat("     CpG/Gene location file file: ",geneLocation.csv.file,"\n")
cat("     Factor covariates: ",covar.fact,"\n")
cat("     Numeric covariates: ",covar.num,"\n")
cat("     Chromosome: ",chr,"\n")
cat("     Output files prefix: ",OutPrefix,"\n")
cat("     Overwrite previous results? ",ifelse(overwrite,"Yes","No"),"\n")
cat("#################################################################################################\n")

if(chr=="all"){
  chr <- seq(1,22,1)
}else{
  chr <- as.numeric(str_split(trimws(chr), pattern = ",", simplify = T)[1,])
}

exp.txt.file = paste0(OutPrefix,".exp.txt")
geneLocation.txt.file = paste0(OutPrefix,".gene.loc.txt")
covariat.file = paste0(OutPrefix,".covariates.txt")

if(file.exists(fam.file)){
  cat("Reading fam file...\n")
  samples <- read.table(file=fam.file,sep=" ",header=F,stringsAsFactors=F)
  rownames(samples) <- samples$V2
}else{
  stop("Unable to find ",fam.file)
}

if(file.exists(eigenvec.file)){
  cat("Reading genotype eigenvectors...\n")
  eigenvec <- read.table(file=eigenvec.file,header=F,row.names = 2)[,c(-1)]
}else{
  stop("Unable to find ",eigenvec.file)
}

if(file.exists(exp.rds.file)){
  cat("Reading expression data...\n")
  exp_all <- as.data.frame(readRDS(exp.rds.file))
}else{
  stop("Unable to find ",exp.rds.file)
}

if(!all(sapply(list(colnames(exp_all),rownames(eigenvec)), FUN = identical, rownames(samples)))){
  warning("Sample names in expression matrix and genotype data are not exactly matched! The intersection will be consider.")
  shared_names <- Reduce(intersect, list(colnames(exp_all),rownames(samples),rownames(eigenvec)))
  if(length(shared_names)==0){
    stop("There is no shared samples between the Methylation/Expression matrix and Genotype data! 
         Check the sample IDs in fam and Expression/Methylation files.")
  }else{
    exp_all <- exp_all[,colnames(exp_all) %in% shared_names]
    samples <- samples[rownames(samples) %in% shared_names,]
    
    exp_all <- exp_all[shared_names]
    samples <- samples[match(shared_names,rownames(samples)),]
    eigenvec <- eigenvec[match(shared_names,rownames(eigenvec)),]
  }
  
}

if(file.exists(geneLocation.csv.file)){
  cat("Reading Gene location data...\n")
  gene_loc_all <- read.csv(file=geneLocation.csv.file,stringsAsFactors = F)[,1:4]
  names(gene_loc_all) <- c("geneid","chr","left","right")
}else{
  stop("Unable to find ",geneLocation.csv.file)
}
if(!all(gene_loc_all$geneid %in% row.names(exp_all))){
  warning("Some Gene/CpG IDs in the Gene/CpG location data are not represented in the expression matrix! Removing them...")
  gene_loc_all = gene_loc_all[gene_loc_all$geneid %in% rownames(exp_all),]
}

exp_all <- cbind.data.frame(geneid = rownames(exp_all),exp_all)


if((!file.exists(exp.txt.file))|(overwrite)){
  cat("Saving Expression data...\n")
  write.table(exp_all,file = exp.txt.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(exp.txt.file," exists\n")
}
if((!file.exists(geneLocation.txt.file))|(overwrite)){
  cat("Saving gene location data...\n")
  write.table(gene_loc_all,file = geneLocation.txt.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(geneLocation.txt.file," exists\n")
}

for (i in chr) {
  cat("************************\n")
  cat("Working on chr ",i,":\n")
  
  genotype.txt.file = paste0(OutPrefix,".chr",i,".genotype.txt")
  genotype.raw.file = paste0(OutPrefix,".chr",i,".raw")
  snp.file = paste0(OutPrefix,".snps",".chr",i,".txt")
  snp.loc.file=paste0(OutPrefix,".snps.loc",".chr",i,".txt")
  
  if((!file.exists(snp.file)) | (!file.exists(snp.loc.file))|(overwrite)){
    if(file.exists(genotype.txt.file)){
      cat(paste("Reading Genotype text data chr",i,"...\n"))
      snps<-fread(genotype.txt.file, data.table = FALSE,nThread=16)
    } else
      if(file.exists(genotype.raw.file)){
        
        cat(paste("Reading Genotype raw data chr",i,"...\n"))
        snps<-fread(genotype.raw.file, data.table = FALSE,nThread=16)
        
        snps <- t(snps)
        snps <- cbind.data.frame(rownames(snps),snps)
        names(snps) <- c("snpid",snps[2,c(-1)])
        snps <- snps[c(-1:-6),]
        write.table(snps,file=genotype.txt.file,sep="\t",row.names = F,col.names = T,quote = F)
      }else{
        stop("Unable to find genotype files for chr ",i,":\n",genotype.txt.file,"\n",genotype.raw.file)
      }
    
    if(!identical(colnames(exp_all),colnames(snps)[-1])){
      index <- match(colnames(exp_all)[-1],colnames(snps)[-1])
      snps <- snps[,c(1,index+1)]
    }
    snps_loc <- str_split(snps$snpid,pattern=':',simplify=T)[,c(1,2)]
    snps_loc <- cbind.data.frame(snps$snpid,snps_loc)
    names(snps_loc) <- c("snpid","chr","pos")
    
    cat(paste("Are the SNP ids matched in in chr",i," data? "))
    cat(ifelse(identical(snps$snpid,snps_loc$snpid),"Yes","NO"),"\n")
    
    cat(paste("Are the sample ids matched in expression data and chr",i," SNP data? "))
    cat(ifelse(identical(colnames(exp_all)[-1], colnames(snps)[-1]),"Yes","NO"),"\n")
    
    cat(paste("Saving Chr",i,"...\n"))
    write.table(snps,file = snp.file,quote = F,col.names = T,row.names = F,sep = '\t')
    write.table(snps_loc,file = snp.loc.file,quote = F,col.names = T,row.names = F,sep = '\t')
    
    
  }else
    cat("Both files already exist:\n",snp.file,"\n",snp.loc.file,"\n")
  
}
cat("******************************************\n")
if((!file.exists(covariat.file))|(overwrite)){
  
  cat("Generating covariates file...\n")
  
  covariates <- cbind.data.frame(eigenvec$V3,eigenvec$V4,eigenvec$V5)
  names(covariates) <- c("PC1","PC2","PC3")
  
  if(file.exists(exp.pheno.file)){
    exp.pheno <- read.csv(exp.pheno.file,row.names=1,stringsAsFactors = F) 
    
    if(!identical(rownames(exp.pheno) , colnames(exp_all)[-1])){
      warning("Sample names in expression matrix and phenotype data are not exactly matched! The intersection will be consider.")
      shared_names <- intersect(colnames(exp_all),rownames(exp.pheno))
      exp.pheno <- exp.pheno[rownames(exp.pheno) %in% shared_names,]
      exp.pheno <- exp.pheno[match(shared_names,rownames(exp.pheno)),]
    }
    if(covar.fact!=""){
      covar.fact = trimws(str_split_1(covar.fact,pattern = ','))
      covar.fact = covar.fact[covar.fact != ""]
      if(length(covar.fact) > 0){
        for (c in covar.fact) {
          exp.pheno[,c] <- as.numeric(as.factor(exp.pheno[,c]))
        }
        covariates <- cbind.data.frame(exp.pheno[,covar.fact],covariates)
      }
    }
    if(covar.num != ""){
      covar.num = trimws(str_split_1(covar.num,pattern = ','))
      covar.num = covar.num[covar.num != ""]
      if(length(covar.num) > 0){
        for (c in covar.num) {
          exp.pheno[,c] <- as.numeric(exp.pheno[,c])
        }
        covariates <- cbind.data.frame(exp.pheno[,covar.num],covariates)
      }
    }
  }else{
    message("No phenotype file is provided!")
  }
  
  rownames(covariates) <- colnames(exp_all)[-1]
  covariates <- t(covariates)
  covariates <- cbind.data.frame(rownames(covariates),covariates)
  names(covariates)[1] <- "id"
  
  cat("Sample ids matched in covariate file with the others? ")
  if(file.exists(exp.pheno.file)){
    cat(ifelse(all(sapply(list(colnames(exp_all)[-1], rownames(exp.pheno),rownames(eigenvec),colnames(covariates)[-1]), FUN = identical, rownames(samples))),"Yes","NO"),"\n")
  } else{
    cat(ifelse(all(sapply(list(colnames(exp_all)[-1],rownames(eigenvec),colnames(covariates)[-1]), FUN = identical, rownames(samples))),"Yes","NO"),"\n")
  }
  
  cat("Saving covaraites file...\n")
  write.table(covariates,file = covariat.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(covariat.file," exist\n")
}

cat("Preparing data is done!\n")
