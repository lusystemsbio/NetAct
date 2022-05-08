############################################################################
# To generate network from Data (Using microarray data)
############################################################################

library(NetAct)
library(Biobase)

data("GSE17708")
# Here GSE17708 is the data from GEO database. 

data("gsdb/hDB")
# This is the gene interaction database. Use mDB for mice. 

compList= c("Early-Middle", "Middle-Late", "Early-Late")
# The list of comparison between different datatypes

eset = GSE17708

DErslt = MultiMicroDegs(eset = eset, compList = compList)
# Find differentially expressed genes

GSEA1 = TF_GSEA(GSDB = hDB, DErslt = DErslt$`Early-Middle`, nperm = 1000, qval=T)
GSEA2 = TF_GSEA(GSDB = hDB, DErslt = DErslt$`Middle-Late`, nperm = 1000, qval=T)
GSEA3 = TF_GSEA(GSDB = hDB, DErslt = DErslt$`Early-Late`, nperm = 1000, qval = T)
save(GSEA1, GSEA2, GSEA3, file = "GSEArslt/GSE17708_GSEA.RData")
# GSEA result for different datatypes

#load("GSEArslt/GSE17708_GSEA.RData")

# Identify the DE transcription factors. If some important tfs are missing, use a higher threshold, for example, change 0.05 to 0.1
tfs = unique(c(rownames(GSEA1)[GSEA1$qvals <0.05], rownames(GSEA2)[GSEA2$qvals < 0.05],
               rownames(GSEA3)[GSEA3$qvals < 0.05]))

# Calucalte the activity of each tf
acts = TF_Activity(tfs, hDB, GSE17708, DErslt$Overall)
acts = acts$all_activities

# Plot the activity and expression of tfs
Combine_heatmap(acts, GSE17708)

exprs = exprs(GSE17708)[rownames(acts),]
save(acts, exprs, file = "Acts/GSE17708_Racipe.RData")

# Construct the network. adjust the network size using miTh. Network will be in network.txt file
tf_links <- TF_Filter(acts, GSDB = hDB, miTh = 0.55, nbins = 8, method = "spearman")

# Tfs in the network
tfsInNetwork <- sort(unique(union(tf_links$from, tf_links$to)))
tfsInNetwork

# Plot the network. Check if the network is broken. 
plot_network_v(tf_links)


############################################################################
# RNAseq example, Use voom to convert RNAseq to microarray type data

# read counts data
data_emt <- read.table("data/GSE17708.RData", header = T)

grouping <- c (rep("C", each = 8), rep("EMT", each = 8))
# convert symbols

t = sapply(data_emt$gene_id, function(x) strsplit(x, "[.]")[[1]], USE.NAMES = FALSE)
data_emt$gene_id = t[1,]
data_emt <- data_emt[!duplicated(data_emt$gene_id),]
counts <- data_emt[,3:18]

geneid = data_emt$gene_id
genes <- select(org.Hs.eg.db, keys=geneid, columns="SYMBOL", keytype="ENSEMBL")
genes = genes[!duplicated(genes$ENSEMBL),]

y <- DGEList(counts=counts,group=grouping, genes = genes)

keep <- rowSums(y$count>10) >= 1
hasnona <- rowSums(is.na(y$genes)) == 0
hasnodup <- !duplicated(y$genes$SYMBOL)
y <- y[keep & hasnona & hasnodup, keep.lib.size = FALSE]
nrow(y$counts)

x <- calcNormFactors(y)
x$sample

lcpm = cpm(x, normalized.lib.sizes=TRUE, log=TRUE)
boxplot(lcpm)
pca = prcomp(t(lcpm), scale = TRUE)
plot(pca$x[,1:2])
text(pca$x[,1:2]+3, labels = grouping, cex = 0.5)
pca$sdev

###### deg with limma
compare_list = c("C-EMT")

library(limma)
ct = factor(grouping)
design <- model.matrix(~0 + ct)
colnames(design) <- levels(ct)
contr.matrix <- makeContrasts(contrasts = compare_list,
                              levels = colnames(design))

v <- voom(x, design, plot = TRUE)
boxplot(v$E)
#edata <- as.data.frame(v$E)
#rownames(edata) = v$genes$SYMBOL
e = v$E
rownames(e) = v$genes$SYMBOL

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
results = decideTests(efit, p.value = 0.05)
summary(results)

fit2 <- treat(efit,lfc=0.1)
volcanoplot(fit2,coef=1)

#### save data for NetAct
rslt = topTable(efit, coef=1, number = Inf, sort.by = "P")
rslt$padj = rslt$adj.P.Val
rownames(rslt) = rslt$SYMBOL
rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
degs = rownames(rslt[rslt$padj< 0.05, ])
DErslt = list(table = rslt, rank_vector = rank_vector, degs = degs)

phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = grouping))
rownames(phenoData) = colnames(e)
neweset = ExpressionSet(assayData = e, phenoData = phenoData)
save(neweset, DErslt, file = "emt_adam.RData")

#### NetAct
library(NetAct)
data("hDB_v2");data("hDB")
load("emt_adam.RData")
#gsearslt_v2 = TF_GSEA(hDB_v2, DErslt, minSize=8, nperm = 10000, qval = T)
#write.csv(gsearslt_v2, "results_hDB_v2.csv")
#gsearslt = TF_GSEA(hDB, DErslt, minSize=8, nperm = 10000, qval = T)
#write.csv(gsearslt, "results_hDB.csv")

gsearslt = read.csv("results_hDB.csv")
tfs = gsearslt$tf[gsearslt$qvals <0.2]

acts_mat = TF_Activity(tfs, hDB, neweset, DErslt)$all_activities

