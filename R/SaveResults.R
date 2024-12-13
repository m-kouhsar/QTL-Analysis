args<-commandArgs(trailingOnly = TRUE)

OutPrefix <- args[1]
cis_pval <- as.numeric(args[2])
trans_pval <- as.numeric(args[3])
dist_ <- as.numeric(args[4])
trans.cross.chr <- ifelse(tolower(trimws(args[5]))=="yes",".CrossChr","")
save.csv.cis <- ifelse(tolower(trimws(args[6]))=="yes",T,F)
save.csv.trans <- ifelse(tolower(trimws(args[7]))=="yes",T,F)
chr <- trimws(args[8])
overwrite <- ifelse(trimws(tolower(args[9]))=="yes" , T ,F)

library(stringr)

if(chr=="all"){
  chr <- seq(1,22,1)
}else{
  chr <- as.numeric(str_split(trimws(chr), pattern = ",", simplify = T)[1,])
}

OutDir <- dirname(OutPrefix)
ResultsDir <- paste0(OutDir , "/QTL.Results")
setwd(ResultsDir)
OutFilePrefix <- basename(OutPrefix)

cis_file = paste0(OutFilePrefix,".eQTL.Cis.Pval.",cis_pval,".Dist.",dist_,".csv")
trans_file = paste0(OutFilePrefix,".eQTL.Trans",trans.cross.chr,".Pval.",trans_pval,".Dist.",dist_,".csv")
merge_file = paste0(OutFilePrefix,".eQTL",trans.cross.chr,".TransPvalue.",trans_pval,".CisPvalue.",cis_pval,".Dist.",dist_,".RData")

cat("Preparing merged QTL results...\n")

if((!file.exists(merge_file))|(overwrite)){
  for (i in chr) {
    file_ <- paste0(OutFilePrefix,".matrixEQTL.chr",i,".RData")
    if(file.exists(file_)){
      cat("Reading QTL file Chr",i,"(",file_,")","...\n")
      load(file_)
      assign(x = paste0("eQTL.chr",i),value = me)
      remove(me)
    }else{
      cat("There is no results for chromosome ",i,"\n")
    }
    
  }
  
  variables=ls()[str_detect(ls(),pattern ="eQTL.chr" )]
  mQTLs <- mget(variables)
  cat("Writing merged rdat file...\n")
  save(mQTLs,file = merge_file)
}else{
  cat("Loading merged rdat file...\n")
  load(merge_file)
}

mqtl.cis <- vector(mode = "list",length = length(mQTLs))
mqtl.trans <- vector(mode = "list",length = length(mQTLs))
for (i in 1:length(mQTLs)) {
  if(save.csv.cis)
    mqtl.cis[[i]] <- mQTLs[[i]]$cis$eqtls[mQTLs[[i]]$cis$eqtls$pvalue	< cis_pval,]
  if(save.csv.trans)
    mqtl.trans[[i]] <- mQTLs[[i]]$trans$eqtls[mQTLs[[i]]$trans$eqtls$pvalue	< trans_pval,]
}

if(save.csv.cis)
  mqtl.cis_out <- do.call(rbind.data.frame,mqtl.cis)
if(save.csv.trans)
  mqtl.trans_out <- do.call(rbind.data.frame,mqtl.trans)

if(save.csv.cis)
  if(ncol(mqtl.cis_out) > 0){
    names(mqtl.cis_out)[2] <- "gene"
    mqtl.cis_out <- mqtl.cis_out[order(mqtl.cis_out$pvalue , decreasing = F),]
    cat("Writing Cis csv file...\n")
    write.csv(mqtl.cis_out,file = cis_file,row.names = F)
  }else{
    cat("There is no Cis results!\n")
  }

if(save.csv.trans)
  if(ncol(mqtl.trans_out) > 0){
    names(mqtl.trans_out)[2] <- "gene"
    mqtl.trans_out <- mqtl.trans_out[order(mqtl.trans_out$pvalue , decreasing = F),]
    cat("Writing Trans csv file...\n")
    write.csv(mqtl.trans_out,file =trans_file,row.names = F)
  }else{
    cat("There is no Trans results!\n")
  }

