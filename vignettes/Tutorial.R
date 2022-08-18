## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 7, fig.height = 6, fig.align = "center")

## -----------------------------------------------------------------------------
library(edgeR)
library(NetAct)
library(sRACIPE)

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")

x <- readDGE(files, columns=c(1,3))
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
x$samples$group <- group
samplenames <- c("LP1", "ML1", "Basal1", "Basal2", "ML2", "LP2", "Basal3", "ML3", "LP3")
colnames(x) <- samplenames

## -----------------------------------------------------------------------------
compList <- c("Basal-LP", "Basal-ML", "LP-ML")
phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = group))
rownames(phenoData) = colnames(x$counts)

## -----------------------------------------------------------------------------
counts <- preprocess_counts(counts = x$counts, groups = group, mouse = TRUE)
DErslt = RNAseqDegs_limma(counts = counts, phenodata = phenoData, 
                          complist = compList, qval = 0.05)
neweset = ExpressionSet(assayData = as.matrix(DErslt$Overall$e), phenoData = phenoData)

## -----------------------------------------------------------------------------
data("mDB")
calc <- FALSE

if (calc) {
  gsearslts <- TF_Selection(GSDB = mDB, DErslt = DErslt, minSize = 5, nperm = 10000,
                            qval = 0.05, compList = compList,
                            nameFile = "gsearslts_tutorial")
} else {
  gsearslts <- readRDS(file = "gsearslts_tutorial.RDS")
}

tfs <- gsearslts$tfs
tfs

## -----------------------------------------------------------------------------
Reselect_TFs(GSEArslt = gsearslts$GSEArslt, qval = 0.01)

## -----------------------------------------------------------------------------
act.me <- TF_Activity(tfs, mDB, neweset, DErslt$Overall)
acts_mat = act.me$all_activities
Activity_heatmap(acts_mat, neweset)

## -----------------------------------------------------------------------------
tf_links = TF_Filter(acts_mat, mDB, miTh = .05, nbins = 8, corMethod = "spearman", DPI = T)
#plot_network(tf_links)

## -----------------------------------------------------------------------------
racipe_results <- sracipeSimulate(circuit = tf_links, numModels = 200, plots = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

