rm(list = ls())
setwd("/Users/koharv/Documents/Work/NetAct/")
load("skd_expression/GSE31912eset.RData")
counts <- exprs(GSE31912eset)
library(NetAct)
###############################################################################
pd <- pData(GSE31912eset)
Genes <- rownames(counts)
DErslt <- list()
padj <- rep(1,length(Genes))
table <- data.frame(padj)
DErslt$table <- table
rownames(DErslt$table) <- rownames(counts)
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
acts_mat = TF_Activity(names(GSDB), GSDB, counts, DErslt, with_weight = F)
saveRDS(acts_mat, "GSE31912_Activity_NetAct.RDS")
GSE31912ActivityNetAct <- acts_mat$all_activities
colnames(GSE31912ActivityNetAct) <- pd$celltype
heatmap.2(as.matrix(GSE31912ActivityNetAct), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "s"))))

kdGene <- pd$celltype
GSE31912ActivityNetAct <- GSE31912ActivityNetAct[which(rownames(GSE31912ActivityNetAct) %in% kdGene),]

heatmap.2(as.matrix(GSE31912ActivityNetAct), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "s"))))
zScoreNetAct <- as.list(rownames(GSE31912ActivityNetAct))
names(zScoreNetAct) <- rownames(GSE31912ActivityNetAct)
#z <- matrix(nrow = length(rownames(tmp)), ncol = 3)
#names(z) <- rownames(tmp)
for(i in 1:dim(GSE31912ActivityNetAct)[1]){
  tf = rownames(GSE31912ActivityNetAct)[i]
  print(tf)
  tmpData <- GSE31912ActivityNetAct[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (tmp[kdNumber] - tmpMean)/tmpSd
  zScoreNetAct[[i]] <- (GSE31912ActivityNetAct[i,kdNumber] - tmpMean)/tmpSd
}
zScoreNetAct
tmp <- unlist((zScoreNetAct))
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
tmp <- data.frame(names(tmp),tmp)
write.table(tmp,file = "zScore_NetAct_GSE31912.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_NetAct_GSE31912.RDS")

####### Viper activity
hDB = hDB_v2
# hDB = hDB1
# 1. Bootstrap + ARANCE + Viper
regulon1 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE31912_mi.txt", counts)
# viper1 = viper(counts, regulon1, verbose = TRUE, minsize = 5, pleiotropy = T)
# saveRDS(viper1, "GSE31912_Activity_Viper5.RDS")
viper1 = viper(counts, regulon1, verbose = TRUE, minsize = 5, pleiotropy = T)
saveRDS(viper1, "GSE31912_Activity_Viper5.RDS")
# viper1 = exprs(viper1)
tfs = rownames(viper1) # as.character(pd$celltype[pd$celltype %in% rownames(viper1)])
Acts1 = viper1[tfs, ]
colnames(Acts1) <- pd$celltype

heatmap.2(as.matrix(Acts1), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x)))))

heatmap.2(as.matrix(tmp), trace = "none")

zScoreViper <- as.list(rownames(Acts1))
names(zScoreViper) <- rownames(Acts1)
#z <- matrix(nrow = length(rownames(tmp)), ncol = 3)
#names(z) <- rownames(tmp)
for(i in 1:dim(Acts1)[1]){
  tf = rownames(Acts1)[i]
  print(tf)
  tmpData <- Acts1[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (tmp[kdNumber] - tmpMean)/tmpSd
  zScoreViper[[i]] <- (Acts1[i,kdNumber] - tmpMean)/tmpSd
}
zScoreViper
tmp <- unlist(zScoreViper)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
tmp <- data.frame(names(tmp),tmp)
write.table(tmp,file = "zScore_ViperAracne_GSE31912.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_ViperAracne_GSE31912.RDS")

# GSDB+VIper
tf <- rownames(GSE31912ActivityNetAct)
tgt <- unlist(sapply(tf, function(x) unlist(GSDB[[x]])))
tf <- unlist(sapply(tf, function(x) rep(x,length(GSDB[[x]])),USE.NAMES = F))
aracneGSDB <-data.frame(tf = tf,target = tgt,mi=rep(1,length(tf)))
write.table(aracneGSDB, file = "skd_expression/outputFolder/Bootstrap_GSE31912_mi_GSDB.txt", row.names = F,quote = F, sep = "\t")
# write gene as the first letter in this file

regulon3 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE31912_mi_GSDB.txt", eset = counts)
viper3 = viper(counts, regulon3, verbose = T, minsize = 5, pleiotropy = T)
Acts3 = viper3
colnames(Acts3) <- pd$celltype

heatmap.2(as.matrix(Acts3), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "spearman"))))

heatmap.2( as.matrix(counts), trace = "none")

zScoreViperGSDB <- as.list(rownames(Acts3))
names(zScoreViperGSDB) <- rownames(Acts3)
#z <- matrix(nrow = length(rownames(Acts1)), ncol = 3)
#names(z) <- rownames(Acts1)
for(i in 1:dim(Acts3)[1]){
  tf = rownames(Acts3)[i]
  print(tf)
  
  tmpData <- Acts3[i,]
  #   tmpData <- Acts1[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScoreViperGSDB[[i]] <- (Acts3[i,kdNumber] - tmpMean)/tmpSd
}
zScoreViperGSDB
tmp <- unlist(zScoreViperGSDB)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
tmp <- data.frame(names(tmp),tmp)
write.table(tmp,file = "zScore_ViperGSDBsameMI_GSE31912.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_ViperGSDBsameMI_GSE31912.RDS")

Genes <- rownames(counts)
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
tf <- names(GSDB)
tf <- tfs
colnames <- seq(1:dim(counts)[2])
counts2 <- counts
colnames(counts2) <- colnames
write.table(counts2, file = "skd_expression/expressionGSE31912.txt", row.names = T,quote = F, sep = "\t")
# write gene as the first letter in this file
write.table(unique(tf), file = "skd_expression/tfListGSE31912.txt", row.names = F,quote = F, sep = "\t", col.names = F)
# Generate mutual information using ARACNE code after manually changing the mutual
# information threshold to 1e-8 to ensure that all interactions in the database
# get some mutual information value. Otherwise interactions with MI less than the
# threshold will be missing. Disable DPI as well

# java -Xmx5G -jar dist/aracne.jar -e test/expressionGSE27869.txt  -o outputFolder --tfs test/tfListGSE27869.txt --pvalue 1E-8 --nodpi --seed 1

aracneGSDB <- read.table("skd_expression/outputFolder/single_GSE31912.txt", stringsAsFactors = F, header = T)
aracneGSDB <- aracneGSDB[which(aracneGSDB[,1] %in% tf),]

GSDBNew <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(GSDBNew) <- c("Regulator", "Target", "MI")

for(i in seq_along(tf)){
  tfTmp <- tf[i]
  GSDBtmp <- aracneGSDB[which(aracneGSDB$Regulator == tfTmp),]
  GSDBtmp <- GSDBtmp[which(GSDBtmp$Target %in% GSDB[[tfTmp]]),]
  
  GSDBNew <- rbind(GSDBNew,GSDBtmp)
}

write.table(GSDBNew, file = "skd_expression/outputFolder/Bootstrap_GSE31912_dataMI_GSDB.txt", row.names = F,quote = F, sep = "\t")


regulon3 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE31912_dataMI_GSDB.txt", eset = counts)
viper3 = viper(counts, regulon3, verbose = T, minsize = 5, pleiotropy = T)
Acts3 = viper3
colnames(Acts3) <- pd$celltype

zScoreViperGSDB <- as.list(rownames(Acts3))
names(zScoreViperGSDB) <- rownames(Acts3)
#z <- matrix(nrow = length(rownames(Acts1)), ncol = 3)
#names(z) <- rownames(Acts1)
for(i in 1:dim(Acts3)[1]){
  tf = rownames(Acts3)[i]
  print(tf)
  
  tmpData <- Acts3[i,]
  #   tmpData <- Acts1[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScoreViperGSDB[[i]] <- (Acts3[i,kdNumber] - tmpMean)/tmpSd
}
zScoreViperGSDB
tmp <- unlist(zScoreViperGSDB)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
tmp <- data.frame(names(tmp),tmp)
write.table(tmp,file = "zScore_ViperGSDBdataMI_GSE31912.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_ViperGSDBdataMI_GSE31912.RDS")


###############################################################################
