args<-commandArgs(trailingOnly = TRUE)

OutPrefix <- args[1]
cis_pval <- as.numeric(args[2])
trans_pval <- as.numeric(args[3])
dist_ <- as.numeric(args[4])
trans.cross.chr <- ifelse(tolower(trimws(args[5]))=="yes",".CrossChr","")
save.csv.cis <- ifelse(tolower(trimws(args[6]))=="yes",T,F)
save.csv.trans <- ifelse(tolower(trimws(args[7]))=="yes",T,F)

library(stringr)
OutDir <- dirname(OutPrefix)
setwd(OutDir)
OutPrefix <- basename(OutPrefix)

cis_file = paste0(OutPrefix,".eQTL.Cis.Pval.",cis_pval,".Dist.",dist_,".csv")
trans_file = paste0(OutPrefix,".eQTL.Trans",trans.cross.chr,".Pval.",trans_pval,".Dist.",dist_,".csv")
merge_file = paste0(OutPrefix,".eQTL",trans.cross.chr,".TransPvalue.",trans_pval,".CisPvalue.",cis_pval,".Dist.",dist_,".RData")

cat("\n")
print("Preparing merged QTL results...")
cat("\n")

if(!file.exists(merge_file)){
  for (i in 1:22) {
    file_ <- paste0(OutPrefix,".matrixEQTL.chr",i,"/",paste0(OutPrefix,".matrixEQTL.chr",i,".RData"))
    if(file.exists(file_)){
      print(paste("Reading QTL file Chr",i))
      load(file_)
      assign(x = paste0("eQTL.chr",i),value = me)
      remove(me)
    }
    
  }
  
  variables=ls()[str_detect(ls(),pattern ="eQTL.chr" )]
  mQTLs <- mget(variables)
  print("Writing merged rdat file...")
  save(mQTLs,file = merge_file)
}else{
  print("Loading merged rdat file...")
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
    print("Writing Cis csv file...")
    write.csv(mqtl.cis_out,file = cis_file,row.names = F)
  }else{
    print("There is no Cis results!")
  }

if(save.csv.trans)
  if(ncol(mqtl.trans_out) > 0){
    names(mqtl.trans_out)[2] <- "gene"
    mqtl.trans_out <- mqtl.trans_out[order(mqtl.trans_out$pvalue , decreasing = F),]
    print("Writing Trans csv file...")
    write.csv(mqtl.trans_out,file =trans_file,row.names = F)
  }else{
    print("There is no Trans results!")
  }

