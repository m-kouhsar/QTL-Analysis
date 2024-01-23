args<-commandArgs(trailingOnly = TRUE)

InDir<-args[1]
FilePrefix <- args[2]
cis_pval <- as.numeric(args[3])
trans_pval <- as.numeric(args[4])
dist_ <- as.numeric(args[5])
trans.cross.chr <- ifelse(tolower(args[6])=="yes",".CrossChr","")

library(stringr)
setwd(InDir)

cis_file = paste0(FilePrefix,".eQTL.Cis.Pval.",cis_pval,".Dist.",dist_,".csv")
trans_file = paste0(FilePrefix,".eQTL.Trans",trans.cross.chr,".Pval.",trans_pval,".Dist.",dist_,".csv")
merge_file = paste0(FilePrefix,".eQTL",trans.cross.chr,".TransPvalue.",trans_pval,".CisPvalue.",cis_pval,".Dist.",dist_,".RData")

#print(trans_file)
#print(merge_file)
cat("\n")
print("Preparing merged QTL results...")
cat("\n")

for (i in 1:22) {
  if(file.exists(paste0(FilePrefix,".matrixEQTL.chr",i,"/",paste0(FilePrefix,".matrixEQTL.chr",i,".RData")))){
    print(paste("Reading QTL file Chr",i))
    load(paste0(FilePrefix,".matrixEQTL.chr",i,"/",paste0(FilePrefix,".matrixEQTL.chr",i,".RData")))
    assign(x = paste0("eQTL.chr",i),value = me)
    remove(me)
  }

}

variables=ls()[str_detect(ls(),pattern ="eQTL.chr" )]
mQTLs <- mget(variables)
print("Writing merged rdat file...")
save(mQTLs,file = merge_file)

sig_mqtls_cis <- vector(mode = "list",length = length(mQTLs))
sig_mqtls_trans <- vector(mode = "list",length = length(mQTLs))
for (i in 1:length(mQTLs)) {
  sig_mqtls_cis[[i]] <- mQTLs[[i]]$cis$eqtls[mQTLs[[i]]$cis$eqtls$pvalue	< cis_pval,]
  sig_mqtls_trans[[i]] <- mQTLs[[i]]$trans$eqtls[mQTLs[[i]]$trans$eqtls$pvalue	< trans_pval,]
}

sig_mqtls_cis_out <- do.call(rbind.data.frame,sig_mqtls_cis)
sig_mqtls_trans_out <- do.call(rbind.data.frame,sig_mqtls_trans)

if(ncol(sig_mqtls_cis_out) > 0){
  names(sig_mqtls_cis_out)[2] <- "gene"
  sig_mqtls_cis_out <- sig_mqtls_cis_out[order(sig_mqtls_cis_out$pvalue , decreasing = F),]
  print("Writing Cis csv file...")
  write.csv(sig_mqtls_cis_out,file = cis_file,row.names = F)
}else{
  print("There is no Cis results!")
}

if(ncol(sig_mqtls_trans_out) > 0){
  names(sig_mqtls_trans_out)[2] <- "gene"
  sig_mqtls_trans_out <- sig_mqtls_trans_out[order(sig_mqtls_trans_out$pvalue , decreasing = F),]
  print("Writing Trans csv file...")
  write.csv(sig_mqtls_trans_out,file =trans_file,row.names = F)
}else{
  print("There is no Trans results!")
}
