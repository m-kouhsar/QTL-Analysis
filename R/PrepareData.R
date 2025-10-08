args<-commandArgs(trailingOnly = TRUE)

suppressMessages(library(stringr))
suppressMessages(library(data.table))

fam.file <- trimws(args[1])
eigenvec.file <- trimws(args[2])
exp.rds.file <- trimws(args[3])
exp.pheno.file <- trimws(args[4])
geneLocation.file <- trimws(args[5])
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
cat("     CpG/Gene location file: ",geneLocation.file,"\n")
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

OutDir=dirname(OutPrefix)
OutFilePrefix=basename(OutPrefix)

PlinkDir=paste0(OutDir,"/QTL.PlinkData")
FormattedDataDir=paste0(OutDir,"/QTL.PreparedInput")

dir.create(FormattedDataDir,recursive = T , showWarnings = F)

exp.txt.file = paste0(FormattedDataDir,"/",OutFilePrefix,".exp.txt")
geneLocation.txt.file = paste0(FormattedDataDir,"/",OutFilePrefix,".gene.loc.txt")
covariat.file = paste0(FormattedDataDir,"/",OutFilePrefix,".covariates.txt")

if(file.exists(fam.file)){
  cat("Reading fam file...\n")
  samples <- fread(fam.file,header=F,stringsAsFactors=F, data.table = F)
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

cols_exp_all <- colnames(exp_all)
rows_eigenvec <- rownames(eigenvec)
rows_samples <- rownames(samples)
if(!identical(cols_exp_all, rows_eigenvec) || !identical(cols_exp_all, rows_samples)){
  
  sorted_exp_all <- sort(cols_exp_all)
  sorted_eigenvec <- sort(rows_eigenvec)
  sorted_samples <- sort(rows_samples)
  if(identical(sorted_exp_all, sorted_eigenvec) && identical(sorted_exp_all, sorted_samples)){
    exp_all <- exp_all[, sorted_exp_all]
    eigenvec <- eigenvec[sorted_exp_all, ]
    samples <- samples[sorted_exp_all, ]
  }else{
    stop("Sample names in expression matrix and genotype data are not matched! 
       Checked IID in fam file and column names in expression matrix.")
  }
}

if(file.exists(geneLocation.file)){
  cat("Reading Gene location data...\n")
  gene_loc_all <- read.csv(file=geneLocation.file,stringsAsFactors = F)
  colnames(gene_loc_all) <- tolower(colnames(gene_loc_all))
  gene_loc_all <- gene_loc_all[,c("id","chr","start","end")]
  names(gene_loc_all) <- c("geneid","chr","left","right")
}else{
  stop("Unable to find ",geneLocation.file)
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
  cat(exp.txt.file," already exists.\n")
}
if((!file.exists(geneLocation.txt.file))|(overwrite)){
  cat("Saving gene location data...\n")
  write.table(gene_loc_all,file = geneLocation.txt.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(geneLocation.txt.file," already exists.\n")
}

if((!file.exists(covariat.file))|(overwrite)){
  
  cat("Generating covariates file...\n")
  
  covariates <- cbind.data.frame(eigenvec$V3,eigenvec$V4,eigenvec$V5)
  names(covariates) <- c("PC1","PC2","PC3")
  
  if(file.exists(exp.pheno.file)){
    exp.pheno <- read.csv(exp.pheno.file,row.names=1,stringsAsFactors = F) 
    
    covar.num = trimws(str_split_1(covar.num,pattern = ','))
    covar.num = covar.num[covar.num != ""]
    covar.fact = trimws(str_split_1(covar.fact,pattern = ','))
    covar.fact = covar.fact[covar.fact != ""]
    covar.all = c(covar.num , covar.fact)
    print(paste(covar.all , collapse = ","))
    if(!all(c(covar.num , covar.fact) %in% colnames(exp.pheno))){
      stop("The following covariates are not in the phenotype file column names:\n",
           paste(setdiff(c(covar.num , covar.fact) , colnames(exp.pheno))))
    }
    
    if(!identical(rownames(exp.pheno) , colnames(exp_all)[-1])){
      warning("Sample names in expression matrix and phenotype data are not exactly matched! The intersection will be consider.")
      shared_names <- intersect(colnames(exp_all),rownames(exp.pheno))
      exp.pheno <- exp.pheno[rownames(exp.pheno) %in% shared_names,]
      exp.pheno <- exp.pheno[match(shared_names,rownames(exp.pheno)),]
    }
    
    if(length(covar.all) > 0){
      for (c in covar.all) {
        if(c %in% covar.fact)
          exp.pheno[,c] <- as.numeric(as.factor(exp.pheno[,c]))
        else if(c %in% covar.num)
          exp.pheno[,c] <- as.numeric(as.factor(exp.pheno[,c]))  
      }
      covariates <- cbind.data.frame(exp.pheno[,covar.all],covariates)
      
    }else{
      warning("No covariates provided! Only three first PCs from the genotype data will be used as covariates in the model.")
    }
    
  }else{
    warning("No phenotype file is provided!Only three first PCs from the genotype data will be used as covariates in the model.")
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
  cat(covariat.file," already exists.\n")
}

for (i in chr) {
  cat("************************\n")
  cat("Working on chr ",i,":\n")
  
  genotype.txt.file = paste0(PlinkDir,"/",OutFilePrefix,".chr",i,".genotype.txt")
  genotype.raw.file = paste0(PlinkDir,"/",OutFilePrefix,".chr",i,".raw")
  snp.file = paste0(FormattedDataDir,"/",OutFilePrefix,".snps",".chr",i,".txt")
  snp.loc.file=paste0(FormattedDataDir,"/",OutFilePrefix,".snps.loc",".chr",i,".txt")
  
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
        stop("Unable to find genotype files for chr ",i,":\n",genotype.raw.file)
      }
    
    if(!identical(colnames(exp_all),colnames(snps)[-1])){
      index <- match(colnames(exp_all)[-1],colnames(snps)[-1])
      snps <- snps[,c(1,index+1)]
    }
    snps_loc <- str_split(snps$snpid,pattern=':',simplify=T)[,c(1,2)]
    snps_loc <- cbind.data.frame(snps$snpid,snps_loc)
    names(snps_loc) <- c("snpid","chr","pos")
    snps_loc$pos <- as.numeric(gsub("[^0-9]", "", snps_loc$pos))
    
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

cat("Preparing data is done!\n")
