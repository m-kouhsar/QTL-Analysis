
#####################################################################################################
#                                                                                                   #
#Test local and distand gene-SNP pairs separately and plot Q-Q plots of local and distant p-values  #
#                                                                                                   #
#####################################################################################################

args<-commandArgs(trailingOnly = TRUE)

InDir<-args[1]
chr <- args[2]
FilePrefix <- args[3]
Dist <- args[4]
Dist <- as.numeric(Dist)
cis.pval <- args[5]
cis.pval <- as.numeric(cis.pval)
trans.pval <- args[6]
trans.pval <- as.numeric(trans.pval)
trans.cross.chr <- ifelse(tolower(args[7])=="yes",T,F)

flag1=FALSE
flag2=FALSE

print("Arguments: ")
print(paste0("  Working Dir: ",InDir))
print(paste0("  FilePrefix: ",FilePrefix))
print(paste0("  CHR: ",chr))
print(paste0("  cis distance: ",Dist))
print(paste0("  cis Pvalue: ",cis.pval))
print(paste0("  trans Pvalue: ",trans.pval))
print(paste0("  Get the trans results in cross chromosom mode? ",ifelse(trans.cross.chr,"Yes","No")))

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

setwd(InDir)

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Covariates file name
# Set to character() for no covariates
if(file.exists(paste0(FilePrefix,".covariates.txt"))){
  covariates_file_name = paste0(FilePrefix,".covariates.txt")
}else{
  covariates_file_name = character()
}


# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# Only associations significant at this level will be saved
pvOutputThreshold_cis = cis.pval;
pvOutputThreshold_tra = trans.pval;

# Distance for local gene-SNP pairs
cisDist = Dist;
  
# Genotype file name
SNP_file_name = paste0(FilePrefix,".snps.chr",chr,".txt")

# Gene expression file name
expression_file_name = paste0(FilePrefix,".exp.txt")

# Only associations significant at this level will be saved

snps_location_file_name =paste0(FilePrefix,".snps.loc.chr",chr,".txt") ;

gene_location_file_name = paste0(FilePrefix,".gene.loc.txt");


# Output file name


## Load genotype data
tryCatch(expr = {
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$LoadFile(SNP_file_name);
},
error=function(e){
  flag1 <<- TRUE
})

## Load gene expression data
tryCatch(expr = {
expression = SlicedData$new();
expression$fileDelimiter = "\t";      # the TAB character
expression$fileOmitCharacters = "NA"; # denote missing values;
expression$fileSkipRows = 1;          # one row of column labels
expression$fileSkipColumns = 1;       # one column of row labels
expression$LoadFile(expression_file_name);
},
error=function(e){
  flag2 <<- TRUE
  })
## Load covariates
tryCatch(expr = {
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
},
error=function(e){
  print("An error occured during reading the covariates file")
})
## Run the analysis
if((!flag1)&(!flag2)){
  dir.create(paste0(FilePrefix,".matrixEQTL.chr",chr))
  
  output_file_name_cis =paste0(FilePrefix,".matrixEQTL.chr",chr,paste0("/",FilePrefix,".matrixEQTL.cis.chr",chr,".out"))
  output_file_name_tra = paste0(FilePrefix,".matrixEQTL.chr",chr,paste0("/",FilePrefix,".matrixEQTL.trans.chr",chr,".out"))
  
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  if(!trans.cross.chr){
    index <- which(genepos$chr==chr)
    genepos <- genepos[index,]
    expression$RowReorder(index)
  }
  me = Matrix_eQTL_main(
    snps = snps,
    gene = expression,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE );
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  ## Results:
  save(me,file=paste0(FilePrefix,".matrixEQTL.chr",chr,paste0("/",FilePrefix,".matrixEQTL.chr",chr,".RData")))
  print(paste("result saved in",paste0(FilePrefix,".matrixEQTL.chr",chr,paste0("/",FilePrefix,".matrixEQTL.chr",chr,".RData"))))
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
}else{
  if(flag1)
    print(paste("Loading genotype data for chr",chr,"was failed (probabily there is no genotype data)"))
  if(flag2)
    print(paste("Loading expression data for chr",chr,"was failed (probabily there is no expression data)"))
}
