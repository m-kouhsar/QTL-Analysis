suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(rlang)))
suppressWarnings(suppressMessages(library(performance)))
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

#################################################################################################

## function to remove outliers iteratively 
remove_outlier <- function(df, column, threshold = 4) {
  
  repeat {
    mean_val <- mean(df[[column]], na.rm = TRUE)
    sd_val <- sd(df[[column]], na.rm = TRUE)
    
    new_df <- df[abs(df[[column]] - mean_val) <= threshold * sd_val, ]
    
    if (nrow(new_df) == nrow(df)) {
      break
    }
    
    df <- new_df
  }
  
  return(df)
}

#################################################################################################
# DNAm_file: A tsv file contains normalized methylation values. Rows must be CpG IDs and columns must be Sample IDs
# DNA.annot_file: CpG annotation file. It must contains the following columns:
#                         id: CpG IDs
#                         chr: Chromosome
#                         start: Start position
#                         end: End position
# RNAe_file: A tsv file contains normalized expression values. Rows must be RNA/Gene IDs and columns must be Sample IDs
# RNA.annot_file: RNA/Gene annotation file. It must contains the following columns:
#                         id: RNA/Gene IDs
#                         chr: Chromosome
#                         start: Start position
#                         end: End position
# Pheno_file: A tsv file contains phenotype information about samples
# covars_fact: Categorical covariates (eg. Sex) for adding to the lm model (must be represented in column names in Pheno_file). 
#              The first categorical variable will be used for coloring the Plots 
# covars_num: Numerical covariates (eg. Age) for adding to the lm model (must be represented in column names in Pheno_file) 
# d.cis: Distance from start and end of each RNA/Gene to search for CpGs
# PCs: Numebr of Principal Components (in both methylation and expression data) need to be added to the lm model
# OutPrefix: Result file prefix

#################################################################################################
args <- commandArgs(T)

DNAm_file <- trimws(args[1])
DNA.annot_file <- trimws(args[2])
RNAe_file <- trimws(args[3])
RNA.annot_file <- trimws(args[4])
Pheno_file <- trimws(args[5])
covars_fact <- args[6]
covars_num <- args[7]
d.cis <- trimws(args[8])
PCs <- trimws(args[9])
OutPrefix <- trimws(args[10])


message("Input arguments:")
message("        DNA methylation file: ",DNAm_file)
message("        DNA methylation annotation file: ",DNA.annot_file)
message("        RNA expression file: ",DNAm_file)
message("        RNA annotation file: ",DNA.annot_file)
message("        Phenotype information file: ",Pheno_file)
message("        Factor covariates (will be added to the lm model): ",covars_fact)
message("        Numeric covariates (will be added to the lm model): ",covars_num)
message("        Maximum distance between CpGs and RNAs: ",d.cis)
message("        Number of PCs to add to the lm model: ",PCs)
message("        Output files prefix: ",OutPrefix)

################################################################################################

message("Reading phenotype information...")
Pheno <- read.csv(Pheno_file , row.names = 1 , stringsAsFactors = F)

message("Reading methylation data...")
DNAm <- suppressWarnings(fread(DNAm_file , stringsAsFactors = F, fill=F,data.table = F))

message("Reading expression data...")
RNAe <- suppressWarnings(fread(RNAe_file , stringsAsFactors = F, fill=F , data.table = F))

rownames(DNAm) <- DNAm[,1]
DNAm <- DNAm[,-1]

rownames(RNAe) <- RNAe[,1]
RNAe <- RNAe[,-1]

if(!all(identical(rownames(Pheno) , colnames(DNAm)),identical(rownames(Pheno) , colnames(RNAe)))){
  
  shared_samples <- intersect(rownames(Pheno) , colnames(DNAm))
  shared_samples <- intersect(shared_samples , colnames(RNAe))
  
  if(length(shared_samples) < 2){
    stop("Rownames in Phenotype file are not matched with Colnames in RNA expression and DNA methylation files!")
  }else{
    warning("Only ", length(shared_samples), " samples are matched in Phenotype, DNA methylation and RNA expression files.\n
            The analysis will be restricted to these matched samples")
    Pheno <- Pheno[shared_samples , ]
    DNAm <- DNAm[,shared_samples]
    RNAe <- RNAe[,shared_samples]
  }
}

message("Reading CpG annotation data...")
RNA.annot <- suppressWarnings(fread(RNA.annot_file , stringsAsFactors = F , fill = F , data.table = F))
RNA.annot <- RNA.annot[!duplicated(RNA.annot$id),]

message("Reading RNA annotation data...")
DNA.annot <- suppressWarnings(fread(DNA.annot_file , stringsAsFactors = F , fill = F , data.table = F))
DNA.annot <- DNA.annot[!duplicated(DNA.annot$id),]

covars_fact <- trimws(str_split_1(covars_fact , pattern = ","))
covars_num <- trimws(str_split_1(covars_num , pattern = ","))

if(!all(c(covars_fact,covars_num) %in% colnames(Pheno))){
  stop("The following variables are not presented in Pheno file:\n",
       paste(setdiff(c(covars_fact,covars_num) , colnames(Pheno)),collapse = ", "))
}


for (var_ in covars_num) {
  Pheno[,var_] <- as.numeric(Pheno[,var_])
}
for (var_ in covars_fact) {
  Pheno[,var_] <- as.factor(Pheno[,var_])
}

d.cis <- as.numeric(d.cis)
PCs <- as.numeric(PCs)

#################################################################################################

if(PCs > 0){
  pca <- prcomp(t(DNAm), rank. = PCs)
  pca <- pca$x
  pca <- as.data.frame(scale(pca))
  colnames(pca) <- paste(colnames(pca) , "dnam",sep = ".")
  index <- match(rownames(Pheno) , rownames(pca))
  Pheno <- cbind.data.frame(Pheno , pca[index,])
  
  pca <- prcomp(t(RNAe), rank. = PCs)
  pca <- pca$x
  pca <- as.data.frame(scale(pca))
  colnames(pca) <- paste(colnames(pca) , "rna",sep = ".")
  index <- match(rownames(Pheno) , rownames(pca))
  Pheno <- cbind.data.frame(Pheno , pca[index,])
}

RNA.annot <- RNA.annot[!is.na(RNA.annot$id),]
DNA.annot <- DNA.annot[!is.na(DNA.annot$id),]

RNA.annot <- RNA.annot[RNA.annot$id %in% rownames(RNAe),]
DNA.annot <- DNA.annot[DNA.annot$id %in% rownames(DNAm),]

names(RNA.annot) <- tolower(names(RNA.annot))
names(DNA.annot) <- tolower(names(DNA.annot))

if(nrow(RNA.annot) == 0){
  stop("There is no overlap between RNA IDs in RNA expression and annotation files!")
}

if(nrow(DNA.annot) == 0){
  stop("There is no overlap between CpG IDs in DNA methylation and annotation files!")
}

RNA.annot$window_start = pmax(1 , RNA.annot$start - d.cis)
RNA.annot$window_end = RNA.annot$end + d.cis

rna_gr <- GRanges(
  seqnames = RNA.annot$chr,
  ranges = IRanges(start = RNA.annot$window_start, end = RNA.annot$window_end),
  RNA = RNA.annot$id
)

cpg_gr <- GRanges(
  seqnames = DNA.annot$chr,
  ranges = IRanges(start = DNA.annot$start, end = DNA.annot$end),
  CpG = DNA.annot$id
)

message("Finding CpGs around the gene locations (distance ",d.cis , " bp)...")
overlaps <- findOverlaps(cpg_gr, rna_gr)
overlaps <- paste(
  mcols(cpg_gr)$CpG[queryHits(overlaps)],
  mcols(rna_gr)$RNA[subjectHits(overlaps)],
  sep = ","
)
overlaps <- unique(overlaps)
eQTMs <- as.data.frame(str_split(overlaps , pattern = ",", simplify = T))
colnames(eQTMs) <- c("CpG" , "RNA")

## Select CpGs and transcripts with the selected distance to each other

message("Total number of tests: ", nrow(eQTMs))

## Prepare dataframe for results

eQTMs$outlier = NA
eQTMs$est = NA
eQTMs$std_error = NA
eQTMs$pval = NA
eQTMs$fdr = NA
eQTMs$bonf = NA
eQTMs$rse = NA
eQTMs$vif = NA
eQTMs$outlier = NA
eQTMs$cor = NA
eQTMs$cor_pval = NA

Plots <- vector(mode = "list" , length = nrow(eQTMs))
names(Plots) <- paste(eQTMs$CpG , eQTMs$RNA, sep = "_")

## For loop to remove the outliers, calculate the linear regressions and create the correlation plots		 
message("linear regressions analysis...")
for (i in 1:nrow(eQTMs)) {
  
  message("Running test ",i ,"/",nrow(eQTMs), ": ",eQTMs$CpG[i] , ", " , eQTMs$RNA[i])
  ## extract cpg and transcript for regression
  dnam <- as.data.frame(t(DNAm[eQTMs$CpG[i],]))
  dnam$id <- row.names(dnam)
  rna_e <- as.data.frame(t(RNAe[eQTMs$RNA[i],]))
  rna_e$id <- row.names(dnam)
  
  Pheno$id <- row.names(Pheno)
  
  
  ## merge DNAm, mRNA and covariate data data for regression and plotting
  plot_data <- merge(dnam, rna_e, by= "id")
  plot_data <- merge(plot_data, Pheno, by = "id")
  
  ## remove outliers
  n_with_outlier <- nrow(plot_data)
  plot_data <- remove_outlier(plot_data, 2) ## remove outliers based on DNAm
  plot_data <- remove_outlier(plot_data, 3) ## remove outliers based on mRNA
  n_without_outlier <- nrow(plot_data)
  
  
  ## z-standatization for better comparability of the different data sets 
  if(!all(plot_data[,3] == 0) & !all(plot_data[,2] == 0)){ ## Z-standadization is not possible if all values are 0
    
    plot_data[,2] <- scale(plot_data[,2])
    plot_data[,3] <- scale(plot_data[,3])
    
    ## linear regression with covariates 
    lm_model <- as.formula(paste0(names(plot_data)[3], "~", names(plot_data)[2],"+",
                                  paste0(c(covars_fact , covars_num) , collapse = "+"),"+",
                                  paste0("PC",c(1:PCs),".dnam", collapse = "+"),"+",
                                  paste0("PC",c(1:PCs),".rna", collapse = "+")))
    lm_fit <- lm(lm_model, data = plot_data)
    
    if(!is.na(summary(lm_fit)$coefficients[2,4])){
      eQTMs$est[i] <- summary(lm_fit)$coefficients[2,1]
      eQTMs$std_error[i] <- summary(lm_fit)$coefficients[2,2]
      eQTMs$pval[i] <- summary(lm_fit)$coefficients[2,4]
      eQTMs$rse[i] <- summary(lm_fit)$sigma ## residual standard error
      eQTMs$vif[i] <- check_collinearity(lm_fit)$VIF[1] ## variance inflation factor
      eQTMs$outlier[i] <- n_with_outlier - n_without_outlier ## number of outliers
      
      R <- cor.test(plot_data[,2], plot_data[,3], method = "pearson") ## calculate Pearson correlation coefficient
      R.cor <- round(R$estimate , digits = 2)
      R.pval <- formatC(R$p.value, format = "e", digits = 2)
      eQTMs$cor[i] <- as.numeric(R$estimate)
      eQTMs$cor_pval[i] <- R$p.value
      
      cor_plot <- ggplot(plot_data, aes_string(x = names(plot_data)[2], y = names(plot_data)[3])) +
        geom_point(aes_string(color = all.vars(lm_model)[3])) +	## Color of the dots indicates whether AD or HC
        geom_smooth(method = "lm", se = T, color = "black")+ ## regression line
        annotate("text", x = min(plot_data[,2])+0.75, y = max(plot_data[,3])+0.25, label = paste("Correlation:", R.cor,"\n","P-value:",R.pval))+ ## pearson correlation coefficient
        theme_minimal() +
        labs(color = "Disease status") +
        scale_color_manual(
          values = c("AD" = "#813513", "Control" = "#32584B"))
      
      Plots[[i]] <- cor_plot
      
      # ggsave(
      #   filename = paste0("ROSMAP_corplot_", as.character(d.cis), "_", ## TODO: change file name
      #                     names(plot_data)[2], "_", names(plot_data)[3], ".png"), 
      #   plot = cor_plot,
      #   width = 6, height = 4)
    }
  }
}

## remove invalid p-values so that an FDR correction can be performed
eQTMs <- eQTMs[!is.na(eQTMs$pval),] 
eQTMs <- eQTMs[eQTMs$pval != 0,]

## perform FDR correction
eQTMs$fdr <- p.adjust(eQTMs$pval, method = "fdr")
eQTMs$bonf <- p.adjust(eQTMs$pval, method = "bonferroni")

names(DNA.annot) <- paste0("CpG.", names(DNA.annot))
names(RNA.annot) <- paste0("RNA.", names(RNA.annot))

eQTMs <- merge(eQTMs, DNA.annot, by.x = "CpG", by.y = "CpG.id")
eQTMs <- merge(eQTMs, RNA.annot, by.x = "RNA", by.y = "RNA.id")

eQTMs <- eQTMs[order(eQTMs$pval , decreasing = F),]
## Save results

Plots = Plots[!sapply(Plots, is.null)]

message("Saving the results...")
dir.create(dirname(OutPrefix), showWarnings = F , recursive = T)
save(eQTMs,Plots, file = paste0(OutPrefix,".eQTM.rdat"))

message("All done!")
