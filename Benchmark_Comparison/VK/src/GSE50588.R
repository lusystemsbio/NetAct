setwd("/Users/koharv/Documents/Work/NetAct/")
library(NetAct)
library(gplots)
library(Hmisc)
library(ggplot2)
library(viper)

eset <- readRDS("GSE50588_processed.RDS")
pd <- pData(eset)
counts <- exprs(eset)
pd <- pData(eset)
Genes <- rownames(counts)
DErslt <- list()
padj <- rep(1,length(Genes))
table <- data.frame(padj)
DErslt$table <- table
rownames(DErslt$table) <- rownames(counts)
data("hDB_v2")
#tmp <- readRDS("macrophage/macrophage_DErslt_4h.RDS")
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
kdGene <- pd$`target gene:ch1`
tfs <- names(GSDB)
eset <- counts
acts_mat = TF_Activity(tfs, GSDB, eset, DErslt, with_weight = TRUE, useDatabaseSign = F, useCorSign = T, if_module = T)
# saveRDS(acts_mat, "GSE50588_Activity_NetAct.RDS")
GSE50588ActivityNetAct <- acts_mat$all_activities
colnames(GSE50588ActivityNetAct) <- pd$title
# heatmap.2(as.matrix(GSE50588ActivityNetAct), trace = "none", hclustfun = function(x)
#  hclust(x, method = 'ward.D2'),
#  distfun = function(x)
#    as.dist((1 - cor(t(x), method = "s"))))

GSE50588ActivityNetAct <- GSE50588ActivityNetAct[which(rownames(GSE50588ActivityNetAct) %in% kdGene),]

# heatmap.2(as.matrix(GSE50588ActivityNetAct), trace = "none", hclustfun = function(x)
#  hclust(x, method = 'ward.D2'),
#  distfun = function(x)
#    as.dist((1 - cor(t(x), method = "s"))))
zScoreNetAct <- as.list(rownames(GSE50588ActivityNetAct))
names(zScoreNetAct) <- rownames(GSE50588ActivityNetAct)
#z <- matrix(nrow = length(rownames(tmp)), ncol = 3)
#names(z) <- rownames(tmp)
for(i in 1:dim(GSE50588ActivityNetAct)[1]){
  tf = rownames(GSE50588ActivityNetAct)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- GSE50588ActivityNetAct[i,ctrl]
 # tmpData <- GSE50588ActivityNetAct[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (tmp[kdNumber] - tmpMean)/tmpSd
  zScoreNetAct[[i]] <- (GSE50588ActivityNetAct[i,kdNumber] - tmpMean)/tmpSd
}
zScoreNetAct
tmp <- unlist(zScoreNetAct)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
write.table(tmp,file = "zScore_NetAct_GSE50588.txt", sep="\t",row.names = T,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_NetAct_GSE50588.RDS")

sum(tmp<0)/length(tmp)
sum(tmp<(-1))/length(tmp)
sum(tmp>0)/length(tmp)
sum(tmp>1)/length(tmp)

saveRDS(zScoreNetAct, "zScore_NetAct_GSE50588.RDS")
################# Expression zScore
tfExprs <- counts[which(rownames(counts) %in% kdGene),]
zScoreExprs <- as.list(rownames(tfExprs))
names(zScoreExprs) <- rownames(tfExprs)
#z <- matrix(nrow = length(rownames(tmp)), ncol = 3)
#names(z) <- rownames(tmp)
for(i in 1:dim(tfExprs)[1]){
  tf = rownames(tfExprs)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- tfExprs[i,ctrl]
  tmpData <- tfExprs[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (tmp[kdNumber] - tmpMean)/tmpSd
  zScoreExprs[[i]] <- (tfExprs[i,kdNumber] - tmpMean)/tmpSd
}
zScoreExprs
tmp <- unlist(zScoreExprs)

sum(tmp<0)/length(tmp)
sum(tmp<(-1))/length(tmp)
sum(tmp>0)/length(tmp)
sum(tmp>1)/length(tmp)
#################
####### Viper activity



# eset = GSE27869eset; pd = pData(eset); edata = exprs(eset)
hDB = hDB_v2
# hDB = hDB1
# 1. Bootstrap + ARANCE + Viper
aracne = read.table("skd_expression/outputFolder/Bootstrap_GSE50588_mi.txt", sep="\t", header = T)[, 1:3]
colnames(aracne) = c("tf", "target", "mi")
regulon1 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE50588_mi.txt", counts)
viper1 = viper(counts, regulon1, verbose = TRUE, minsize = 5, pleiotropy = T )
# viper1 = exprs(viper1)
tfs = rownames(viper1) # as.character(pd$celltype[pd$celltype %in% rownames(viper1)])
Acts1 = viper1[tfs, ]
colnames(Acts1) <- pd$title

heatmap.2(as.matrix(Acts1), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x)))))

heatmap.2(as.matrix(tmp), trace = "none")

zScoreViper <- as.list(rownames(Acts1))
names(zScoreViper) <- rownames(Acts1)

for(i in 1:dim(Acts1)[1]){
  tf = rownames(Acts1)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- Acts1[i,ctrl]
#  tmpData <- Acts1[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScoreViper[[i]] <- (Acts1[i,kdNumber] - tmpMean)/tmpSd
}
zScoreViper <- readRDS("zScore_AracneViperGSE50588.RDS")
tmp <- unlist(zScoreViper)
sum(tmp<0)/length(tmp)
sum(tmp<(-1))/length(tmp)
sum(tmp>0)/length(tmp)
sum(tmp>1)/length(tmp)
write.table(tmp,file = "zScore_ViperAracne_GSE50588.txt",sep = "\t", row.names = T,col.names = F, quote = F)

saveRDS(zScoreViper, "zScore_AracneViperGSE50588.RDS")

# GSDB+VIper
tf <- rownames(GSE50588ActivityNetAct)
tgt <- unlist(sapply(tf, function(x) unlist(GSDB[[x]])))
tf <- unlist(sapply(tf, function(x) rep(x,length(GSDB[[x]])),USE.NAMES = F))
aracneGSDB <-data.frame(tf = tf,target = tgt,mi=rep(1,length(tf)))
write.table(aracneGSDB, file = "skd_expression/outputFolder/Bootstrap_GSE50588_mi_GSDB.txt", row.names = F,quote = F, sep = "\t")

regulon3 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE50588_mi_GSDB.txt", eset = counts)
viper3 = viper(counts, regulon3, verbose = T, minsize = 5, pleiotropy = T )
Acts3 = viper3
colnames(Acts3) <- pData(gset)$title

heatmap.2(as.matrix(Acts3), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "spearman"))))

heatmap.2( as.matrix(counts), trace = "none")

zScore <- as.list(rownames(Acts3))
names(zScore) <- rownames(Acts3)
#z <- matrix(nrow = length(rownames(Acts1)), ncol = 3)
#names(z) <- rownames(Acts1)
for(i in 1:dim(Acts3)[1]){
  tf = rownames(Acts3)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- Acts3[i,ctrl]
 # tmpData <- Acts[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScore[[i]] <- (Acts3[i,kdNumber] - tmpMean)/tmpSd
}
# zScore <- readRDS( "zScore_ViperGSDBGSE50588.RDS")
tmp <- unlist(zScore)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
write.table(tmp,file = "zScore_NetAct_GSE50588.txt", sep = "\t",row.names = T,col.names = F, quote = F)


saveRDS(zScore, "zScore_ViperGSDBGSE50588.RDS")
###################################################
Genes <- rownames(counts)
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
tf <- names(GSDB)
counts <- counts[-grep(",",rownames(counts)),]
colnames <- seq(1:dim(counts)[2])
counts2 <- counts
colnames(counts2) <- colnames
write.table(counts2, file = "skd_expression/expressionGSE50588.txt", row.names = T,quote = F, sep = "\t")
# write gene as the first letter in this file
write.table(unique(tf), file = "skd_expression/tfListGSE50588.txt", row.names = F,quote = F, sep = "\t", col.names = F)
# Generate mutual information using ARACNE code after manually changing the mutual
# information threshold to 1e-8 to ensure that all interactions in the database
# get some mutual information value. Otherwise interactions with MI less than the
# threshold will be missing. Disable DPI as well


aracneGSDB <- read.table("skd_expression/outputFolder/single_GSE50588.txt", stringsAsFactors = F, header = T)
aracneGSDB <- aracneGSDB[which(aracneGSDB[,1] %in% tf),]

GSDBNew <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(GSDBNew) <- c("Regulator", "Target", "MI")

for(i in seq_along(tf)){
  tfTmp <- tf[i]
  GSDBtmp <- aracneGSDB[which(aracneGSDB$Regulator == tfTmp),]
  GSDBtmp <- GSDBtmp[which(GSDBtmp$Target %in% GSDB[[tfTmp]]),]
  
  GSDBNew <- rbind(GSDBNew,GSDBtmp)
}

write.table(GSDBNew, file = "skd_expression/outputFolder/Bootstrap_GSE50588_dataMI_GSDB.txt", row.names = F,quote = F, sep = "\t")

regulon3 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE50588_dataMI_GSDB.txt", eset = counts)
viper3 = viper(counts, regulon3, verbose = T, minsize = 5, pleiotropy = T )
Acts3 = viper3
colnames(Acts3) <- kdGene

zScore <- as.list(rownames(Acts3))
names(zScore) <- rownames(Acts3)
#z <- matrix(nrow = length(rownames(Acts1)), ncol = 3)
#names(z) <- rownames(Acts1)
for(i in 1:dim(Acts3)[1]){
  tf = rownames(Acts3)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- Acts3[i,ctrl]
  # tmpData <- Acts[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScore[[i]] <- (Acts3[i,kdNumber] - tmpMean)/tmpSd
}
# zScore <- readRDS( "zScore_ViperGSDBGSE50588.RDS")
tmp <- unlist(zScore)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
write.table(tmp,file = "zScore_ViperGSDBdataMI_GSE50588.txt", sep = "\t",row.names = T,col.names = F, quote = F)


saveRDS(zScore, "zScore_ViperGSDBdataMIGSE50588.RDS")
####################
Acts4 = read.table("skd_expression/NcaActGSE50588.csv", header = T, sep = ",")
rownames(Acts4) <- Acts4$X 
Acts4$X <- NULL
Acts4 <- t(Acts4)

heatmap.2(as.matrix(Acts4), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "spearman"))))

heatmap.2( as.matrix(counts), trace = "none")
Acts4 <- Acts4[kdGene,]
zScore <- as.list(rownames(Acts4))
names(zScore) <- rownames(Acts4)
#z <- matrix(nrow = length(rownames(Acts1)), ncol = 3)
#names(z) <- rownames(Acts1)
for(i in 1:dim(Acts4)[1]){
  tf = rownames(Acts4)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- Acts4[i,ctrl]
  #   tmpData <- Acts1[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScore[[i]] <- (Acts4[i,kdNumber] - tmpMean)/tmpSd
}
zScore
tmp <- unlist(zScore)
sum(tmp>0)/length(tmp)
sum(tmp<0)/length(tmp)
sum(tmp>1)/length(tmp)
sum(tmp<(-1))/length(tmp)


saveRDS(zScore, "zScore_ncaGSE50588.RDS")

####################
zScoreNetAct <- readRDS("zScore_NetAct.RDS")
zScoreViperAracne <- readRDS("zScore_AracneViper.RDS") 
zScoreViperGSDB <- readRDS("zScore_ViperGSDB.RDS")
tfList <- unique(kdGene)
tfList <- tfList[-6]
dataGSE50588 <- list(geneExpression = counts, tfList = tfList)
saveRDS(dataGSE50588, file = "dataGSE50588.RDS")
countNeg <- function(zScore){
  count <- lapply(zScore, function(x){length(which((x< (1)) & (x > (-1)) ))})
  count$zzz <- sum(unlist(count))
  return(count)
}
# Netact, ViperAracne, ViperGSDB
# less than -1: 27/96, 41/96, 19/39
# greater than 1: 27, 36, 15
# between -1 and 1: 42, 19, 5
countNeg(zScoreNetAct)
countNeg(zScoreViperAracne)
countNeg(zScoreViperGSDB)

tfList <- tfList[tfList %in% rownames(counts)]
TfExpr <- counts[tfList,]
zScore <- as.list(rownames(TfExpr))
names(zScore) <- rownames(TfExpr)

for(i in 1:dim(TfExpr)[1]){
  tf = rownames(TfExpr)[i]
  print(tf)
  ctrl <- which(kdGene == "NS")
  tmpData <- TfExpr[i,ctrl]
  #  tmpData <- TfExpr[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (Acts1[kdNumber] - tmpMean)/tmpSd
  zScore[[i]] <- (TfExpr[i,kdNumber] - tmpMean)/tmpSd
}
zScore
###############################################################################
library(GEOquery)
gset <- getGEO("GSE50588", GSEMatrix =TRUE, getGPL=FALSE)
gset <- gset[[1]]
pd <- pData(gset)
data(hDB_v2)
counts <- exprs(gset)
dim(counts)
Genes <- rownames(counts)
# source("https://bioconductor.org/biocLite.R")
# biocLite("illuminaHumanv4.db")

# load library
library(illuminaHumanv4.db)

tmp <- select(illuminaHumanv4.db, 
              keys = Genes, 
              columns=c("SYMBOL"), 
              keytype="PROBEID")


Genes2 <- sapply(Genes, function(x) {tmp$SYMBOL[which(tmp$PROBEID == x)]})
selected <- !duplicated(Genes2)
selected2 <- !is.na(Genes2)
selected <- selected & selected2
counts <- counts[selected,]
dim(counts)
Genes <- Genes2[selected]

rownames(counts) <- Genes
phenoData = new("AnnotatedDataFrame", data = pd)

neweset <- ExpressionSet(assayData = counts, phenoData = phenoData)

saveRDS(neweset, "GSE50588_processed.RDS")
