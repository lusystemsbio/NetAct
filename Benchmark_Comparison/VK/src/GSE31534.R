
rm(list = ls())
setwd("/Users/koharv/Documents/Work/NetAct/")
load("skd_expression/GSE31534eset.RData")
counts <- exprs(GSE31534eset)
library(NetAct)
###############################################################################
pd <- pData(GSE31534eset)
Genes <- rownames(counts)
DErslt <- list()
padj <- rep(1,length(Genes))
table <- data.frame(padj)
DErslt$table <- table
rownames(DErslt$table) <- rownames(counts)
#tmp <- readRDS("macrophage/macrophage_DErslt_4h.RDS")
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
acts_mat = TF_Activity(names(GSDB), GSDB, counts, DErslt, with_weight = F)
GSE31534ActivityNetAct <- acts_mat$all_activities
saveRDS(acts_mat, "GSE31534_Activity_NetAct.RDS")

colnames(GSE31534ActivityNetAct) <- pd$celltype
heatmap.2(as.matrix(GSE31534ActivityNetAct), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "s"))))

kdGene <- pd$celltype
GSE31534ActivityNetAct <- GSE31534ActivityNetAct[which(rownames(GSE31534ActivityNetAct) %in% kdGene),]

heatmap.2(as.matrix(GSE31534ActivityNetAct), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x), method = "s"))))
zScoreNetAct <- as.list(rownames(GSE31534ActivityNetAct))
names(zScoreNetAct) <- rownames(GSE31534ActivityNetAct)
#z <- matrix(nrow = length(rownames(tmp)), ncol = 3)
#names(z) <- rownames(tmp)
for(i in 1:dim(GSE31534ActivityNetAct)[1]){
  tf = rownames(GSE31534ActivityNetAct)[i]
  print(tf)
  tmpData <- GSE31534ActivityNetAct[i,]
  tmpMean <- mean(tmpData)
  tmpSd <- sd(tmpData)
  kdNumber <- which(kdGene == tf)
  
  #  z[i,] <- (tmp[kdNumber] - tmpMean)/tmpSd
  zScoreNetAct[[i]] <- (GSE31534ActivityNetAct[i,kdNumber] - tmpMean)/tmpSd
}
zScoreNetAct
tmp <- unlist(zScoreNetAct)
sum(tmp<0)*100/length(tmp)
sum(tmp< (-1))*100/length(tmp)
tmp <- data.frame(names(tmp),tmp)
write.table(tmp,file = "zScore_NetAct_GSE31534.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_NetAct_GSE31534.RDS")

####### Viper activity
hDB = hDB_v2
# hDB = hDB1
# 1. Bootstrap + ARANCE + Viper
regulon1 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE31534_mi.txt", counts)
viper1 = viper(counts, regulon1, verbose = TRUE, minsize = 5, pleiotropy = T)
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
write.table(tmp,file = "zScore_ViperAracne_GSE31534.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_ViperAracne_GSE31534.RDS")


# GSDB+VIper
tf <- rownames(GSE31534ActivityNetAct)
tgt <- unlist(sapply(tf, function(x) unlist(GSDB[[x]])))
tf <- unlist(sapply(tf, function(x) rep(x,length(GSDB[[x]])),USE.NAMES = F))
aracneGSDB <-data.frame(tf = tf,target = tgt,mi=rep(1,length(tf)))
write.table(aracneGSDB, file = "skd_expression/outputFolder/Bootstrap_GSE31534_mi_GSDB.txt", row.names = F,quote = F, sep = "\t")

regulon3 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE31534_mi_GSDB.txt", eset = counts)
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
write.table(tmp,file = "zScore_ViperGSDBsameMI_GSE31534.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_ViperGSDBsameMI_GSE31534.RDS")

###############################################################################
Genes <- rownames(counts)
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
tf <- names(GSDB)
tf <- tfs
colnames <- seq(1:dim(counts)[2])
counts2 <- counts
colnames(counts2) <- colnames
write.table(counts2, file = "skd_expression/expressionGSE31534.txt", row.names = T,quote = F, sep = "\t")
# write gene as the first letter in this file
write.table(unique(tf), file = "skd_expression/tfListGSE31534.txt", row.names = F,quote = F, sep = "\t", col.names = F)
# Generate mutual information using ARACNE code after manually changing the mutual
# information threshold to 1e-8 to ensure that all interactions in the database
# get some mutual information value. Otherwise interactions with MI less than the
# threshold will be missing. Disable DPI as well

# java -Xmx5G -jar dist/aracne.jar -e test/expressionGSE27869.txt  -o outputFolder --tfs test/tfListGSE27869.txt --pvalue 1E-8 --nodpi --seed 1

aracneGSDB <- read.table("skd_expression/outputFolder/single_GSE31534.txt", stringsAsFactors = F, header = T)
aracneGSDB <- aracneGSDB[which(aracneGSDB[,1] %in% tf),]

GSDBNew <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(GSDBNew) <- c("Regulator", "Target", "MI")

for(i in seq_along(tf)){
  tfTmp <- tf[i]
  GSDBtmp <- aracneGSDB[which(aracneGSDB$Regulator == tfTmp),]
  GSDBtmp <- GSDBtmp[which(GSDBtmp$Target %in% GSDB[[tfTmp]]),]
  
  GSDBNew <- rbind(GSDBNew,GSDBtmp)
}

write.table(GSDBNew, file = "skd_expression/outputFolder/Bootstrap_GSE31534_dataMI_GSDB.txt", row.names = F,quote = F, sep = "\t")


regulon3 = aracne2regulon("skd_expression/outputFolder/Bootstrap_GSE31534_dataMI_GSDB.txt", eset = counts)
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
write.table(tmp,file = "zScore_ViperGSDBdataMI_GSE31534.txt",sep = "\t", row.names = F,col.names = F, quote = F)
saveRDS(zScoreNetAct, "zScore_ViperGSDBdataMI_GSE31534.RDS")



###############################################################################

rm(list = ls())
counts <- read.table("skd_expression/GSE79586_rawcounts_WALSH.txt", header = TRUE, row.names = TRUE)
library(NetAct)
###############################################################################
setwd("~/Documents/Work/NetAct/GSE31534/")
countsData <- read.table(file = "GSE31534_series_matrix.txt", header = FALSE, skip = 63)


library(Biobase)
library(GEOquery)
library(gplots)
# load series and platform data from GEO

gset <- getGEO("GSE31534", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
tmp <- pData(gset$GSE31534_series_matrix.txt.gz)
# set parameters and draw the plot

dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE31534", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)


library(Hmisc)

setwd("/Users/koharv/Documents/Work/NetAct/")
load("skd_expression/GSE27869eset.RData")
data(hDB_v2)
counts <- exprs(GSE27869eset)
dim(counts)
Genes <- rownames(counts)
DErslt <- list()
padj <- rep(1,length(Genes))
table <- data.frame(padj)
DErslt$table <- table
rownames(DErslt$table) <- rownames(counts)
#tmp <- readRDS("macrophage/macrophage_DErslt_4h.RDS")
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
acts_mat = TF_Activity(names(GSDB), GSDB, GSE27869eset, DErslt, with_weight = TRUE)
tmp <- acts_mat$all_activities
colnames(tmp) <- pData(GSE27869eset)$celltype
heatmap.2(as.matrix(Acts4), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x)))))

gsearslt_1 = TF_GSEA(hDB_v2,DErslt, minSize=5, nperm = 500, qval = T, with_we)

pca1 <- prcomp(t((counts)), scale. = TRUE)
dim(pca1$x)
heatmap.2(as.matrix(t(pca1$x[1:51,1:51])), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x)))))



library(Hmisc)

setwd("/Users/koharv/Documents/Work/NetAct/")
load("skd_expression/GSE27869eset.RData")
data(hDB_v2)
counts <- exprs(GSE27869eset)
dim(counts)
Genes <- rownames(counts)
DErslt <- list()
padj <- rep(1,length(Genes))
table <- data.frame(padj)
DErslt$table <- table
rownames(DErslt$table) <- rownames(counts)
#tmp <- readRDS("macrophage/macrophage_DErslt_4h.RDS")
GSDB <- hDB_v2
GSDB <- lapply(GSDB, function(x){if(length(intersect(x,Genes))>4) return(intersect(x,Genes)) else return()})
GSDB <- GSDB[!unlist(lapply(GSDB, is.null))]
acts_mat = TF_Activity(names(GSDB), GSDB, GSE31534eset, DErslt, with_weight = TRUE)
tmp <- acts_mat$all_activities
colnames(tmp) <- pData(GSE31534eset)$celltype
heatmap.2(as.matrix(Acts4), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x)))))

gsearslt_1 = TF_GSEA(hDB_v2,DErslt, minSize=5, nperm = 500, qval = T, with_we)

pca1 <- prcomp(t((counts)), scale. = TRUE)
dim(pca1$x)
pd <- pData(GSE27869eset)
heatmap.2(as.matrix(t(pca1$x[1:51,1:51])), trace = "none", hclustfun = function(x)
  hclust(x, method = 'ward.D2'),
  distfun = function(x)
    as.dist((1 - cor(t(x)))))

###############################################################################