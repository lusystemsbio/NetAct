###########################################################
########### Comparison Between Viper and ARACNE ########### 
###########################################################
rm(list = ls())
library(R.utils)
source("~/Documents/Github/NetAct/R/")
calMI = function(X, Y){
    require(infotheo)
    X1 <- discretize(X, disc="equalfreq")
    Y1 <- discretize(Y, disc="equalfreq")
    mi <- mutinformation(X1,Y1,method="mm")
    return(mi)
}
load("~/Dropbox/NETACT/supplementary/allCombinedDBs.RData")

source("~/Box Sync/Mingyang/logs/Viper_Internal.R")
library(ggplot2)
library(viper)
geoID = "GSE31534eset.RData"
pwy = paste0("~/Box Sync/Mingyang/logs/", geoID)
load(pwy)
eset = GSE31534eset; pd = pData(eset); edata = exprs(eset)
load("~/Box Sync/Mingyang/logs/hgs_Rcistarget_500bp_direct_500_fixed.RData") # load the database
# hDB = test
hDB = hDB1
# 1. Bootstrap + ARANCE + Viper
aracne = read.table("~/Box Sync/Mingyang/ARACNe-AP-master/outputFolder/Bootstrap_GSE31534_mi.txt", sep="\t", header = T)[, 1:3]
colnames(aracne) = c("tf", "target", "mi")
regulon1 = aracne2regulon(aracne, eset)
viper1 = viper(eset, regulon1, verbose = F)
viper1 = exprs(viper1)
tfs = as.character(pd$celltype[pd$celltype %in% rownames(viper1)])
Acts1 = viper1[tfs, ]

# 2. Single Bootstrap + ARACNE + Viper 
aracne = read.table("~/Box Sync/Mingyang/ARACNe-AP-master/outputFolder/Single_GSE31534_mi.txt", sep="\t", header = T)[, 1:3]
colnames(aracne) = c("tf", "target", "mi")
regulon2 = aracne2regulon(aracne, eset)
viper2 = viper(eset, regulon2, verbose = F)
viper2 = exprs(viper2)
tfs = as.character(pd$celltype[pd$celltype %in% rownames(viper2)])
Acts2 = viper2[tfs, ]

# 3. GSDB + Viper
# tt = allNet(hDB)
# tt = tt[ rowSums(apply(tt, 2, function(x) x %in% rownames(eset))) ==2, ]
# miVec = apply(tt, 1, function(x){
#     vec1 = as.vector(exprs(eset)[x[1],])
#     vec2 = as.vector(exprs(eset)[x[2],])
#     return(calMI(vec1, vec2))
#     
# })
# aracne = data.frame(tt, mi = miVec); colnames(aracne) = c("tf", "target", "mi")
load("~/Box Sync/Mingyang/Supplementary/aracne_31534.RData")
regulon3 = aracne2regulon(aracne, eset)
viper3 = viper(eset, regulon3, verbose = F)
viper3 = exprs(viper3)
tfs = as.character(pd$celltype[pd$celltype %in% rownames(viper3)])
Acts3 = viper3[tfs, ]

# 4. GSDB + NETACT 
table = data.frame(padj = rep(1, nrow(eset)))
rownames(table) = rownames(eset)
DErslt = list(table = table)
tmptfs = unique(as.character(pd$celltype[pd$celltype %in% names(hDB)]))
lens = sapply(hDB[tmptfs], function(x) return(length(intersect(x, rownames(eset)))))
tmptfs = names(lens[lens > 4])
netActs = TF_Activity(tmptfs, hDB, eset, DErslt, with_weight = F)
Acts4 = netActs$all_activities

# 5. NCA 
NCA = read.csv("~/Box Sync/Mingyang/Supplementary/nca/rcis.GSE31534.nca.act.csv")
nca = NCA[,2:ncol(NCA)]; rownames(NCA) = NCA$X; nca = t(nca)
colnames(nca) = colnames(Acts4)
Acts5 = nca

Acts = list(Acts1 = Acts1, Acts2 = Acts2, Acts3 = Acts3, Acts4 = Acts4, Acts5 = Acts5)
# save(Acts, file = "~/Box Sync/Mingyang/Supplementary/CompareActs/Acts_Rcis_GSE31534.RData")



# load("~/Box Sync/Mingyang/Supplementary/CompareActs/Acts_Rcis_GSE31534.RData")
geoID = "GSE31534eset.RData"
pwy = paste0("~/Box Sync/Mingyang/logs/", geoID)
load(pwy)
eset = GSE31534eset; pd = pData(eset); edata = exprs(eset)

for (i in 1:length(Acts)){
    tmpActs = Acts[[i]]
    tfs = rownames(tmpActs)
    tfs = intersect(intersect(rownames(tmpActs), pd$celltype), rownames(edata))
    tmpRslt = sapply(tfs, function(tf){
        id = which(pd$celltype == tf)
        exprsVec = as.numeric(edata[tf, ])
        exprsZ = mean(exprsVec[id])- mean(exprsVec)
        actsVec = tmpActs[tf, ]
        tmpMean = mean(actsVec)
        tmpSd = sd(actsVec)
        actsZ = (mean(actsVec[id])- tmpMean)/ tmpSd
        if (i ==4){
            if (actsZ > 0 & exprsZ < 0){
                actsVec = -actsVec
            }
            return((mean(actsVec[id])- mean(actsVec))/sd(actsVec)) # -actsZ
        }else{
            return(actsZ)
        }
    })
    cat("Average Z: ", round(mean(tmpRslt), 2), "\n")
    n1 = sum(tmpRslt< -1); n2 = length(tmpRslt)
    cat(n1, n2, n1/n2, "\n")
}

gs = Reduce(intersect, lapply(Acts, function(x) rownames(x)))
cat(gs, sep = "|")


for (i in 1:length(Acts)){
    tmpActs = Acts[[i]][gs, ]
    tfs = rownames(tmpActs)
    tfs = intersect(rownames(tmpActs), rownames(edata))
    tmpRslt = sapply(tfs, function(tf){
        id = which(pd$celltype == tf)
        exprsVec = as.numeric(edata[tf, ])
        exprsZ = mean(exprsVec[id])- mean(exprsVec)
        actsVec = tmpActs[tf, ]
        tmpMean = mean(actsVec)
        tmpSd = sd(actsVec)
        actsZ = (mean(actsVec[id])- tmpMean)/ tmpSd
        if (i ==4){
            if (actsZ > 0 & exprsZ < 0){
                actsVec = -actsVec
            }
            return((mean(actsVec[id])- mean(actsVec))/sd(actsVec)) # -actsZ
        }else{
            return(actsZ)
        }
        # return(actsZ)
    })
    cat("Average Z: ", round(mean(tmpRslt), 2), "\n")
    n1 = sum(tmpRslt< -1); n2 = length(tmpRslt)
    cat(n1, n2,   "\n")
}




##### Draw the plot #####

# actdata = Acts$Acts4
# tfs = rownames(actdata)
# 
# actdata= t(sapply(tfs, function(tf){
#     id = which(pd$celltype == tf)
#     actsVec = actdata[tf, ]
#     tmpMean = mean(actsVec)
#     if (actsVec[id] > tmpMean){
#         actsVec = -actsVec
#     }
#     return(actsVec)
# }))
# tmpacts = unname(sapply(tfs, function(tf){
#     id = which(pd$celltype == tf)
#     acts = actdata[tf, id]
#     return(acts)}))
# 
# tmpexprs = unname(sapply(tfs, function(tf){
#     id = which(pd$celltype == tf)
#     exprs = edata[tf, id]
#     return(exprs)}))
# 
# pointdata = data.frame(tf = rownames(actdata), acts =  tmpacts, exprs = tmpexprs)
# actdata = data.frame(acts = as.vector(actdata), tf = rep(rownames(actdata), ncol(actdata)))
# ggplot(actdata, aes(x = tf,  y = acts)) + 
#     geom_violin() + geom_boxplot(width = 0.1, fill = "gray90") + 
#     geom_point(data = pointdata, shape = 8, color="red", size = 3,
#                mapping = aes(x = tf, y = acts)) + 
#     theme_classic(base_size = 14) + xlab("")+ ylab("Activities") +
#     ggtitle("RcisTarget DB + NETACT") +
#     theme(legend.position="none",
#           plot.title = element_text(size=14,hjust = 0.5,face="bold"),
#           axis.text.x = element_text(angle = 45, hjust = 1, size=9))
# 






