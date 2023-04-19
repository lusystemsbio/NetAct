##################################################################
################### DE gene analysis functions ###################
##################################################################

#' @title RNA-seq data pre processing
#' @description NetAct uses edgeR to load the count data and the group information for experimental conditions, 
#'              It also coverts gene symbols and remove duplicates.
#' @param counts raw count matrix
#' @param groups group information for experimental conditions
#' @param mouse use mouse genome or not (default: FALSE)
#' @return x$counts: processed count matrix
#' @export
#' @import edgeR
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
Preprocess_counts <- function(counts, groups, mouse = FALSE) {
#  require(edgeR)
#  Create DGEList and Filter
  x <- DGEList(counts, group = groups)
 
#  Gene Annotations and Filter Duplicate Symbols
  geneid <- rownames(counts)
  
  if(grepl("ENS", geneid[1])) {
    keytype = "ENSEMBL"
  } else {
    keytype = "ENTREZID"
  }
  
  if (mouse == TRUE) {
    genes <- AnnotationDbi::select(org.Mm.eg.db, keys = geneid, columns = "SYMBOL", keytype = keytype)
  } else {
    genes <- AnnotationDbi::select(org.Hs.eg.db, keys = geneid, columns = "SYMBOL", keytype = keytype)
  }
  
  dup_id = which(duplicated(genes[keytype]))
  if(length(dup_id) != 0){
    genes = genes[-dup_id,]
  }
  
  filter_na <- which(is.na(genes$SYMBOL))
  filter_dup <- which(duplicated(genes$SYMBOL))
  filter <- unique(filter_na, filter_dup)
  if(length(filter) != 0){
    x$counts <- x$counts[-filter, ]
    rownames(x$counts) <- genes$SYMBOL[-filter]
  }
  
#Filter Low Expression
  keep.exprs <- filterByExpr(x, group=groups)
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  
  return(x$counts)
}

#' @title Helper Function For DEG Analysis of microArray Data (for a single comparison)
#' @param eset Processed gene expression data in the ExpressionSet format
#'              batch & experimental conditions are provided in pData. 
#' @param qval q-value cutoff for DEG analysis (default: 0.05)
#' @return DEG result in the format of a list containing:
#'         table: table of DEG results.
#'         rank_vector: a vector of t-statistics for every gene.
#'         degs: a vector of gene symbols for DEGs.
#' @import edgeR
#' @import limma
DEG_Analysis_Micro <- function(eset, qval = 0.05) {
  data = exprs(eset)
  phenodata = pData(eset)
  celltype = phenodata$celltype
  
  if (is.null(phenodata$batch)){
    fit = lmFit(data, model.matrix(~ celltype))
  }else{
    batch = phenodata$batch
    fit = lmFit(data, model.matrix(~ batch + celltype))
  }
  efit = eBayes(fit)
  rslt = topTable(efit, coef=2, number = Inf, sort.by = "P")
  rslt$padj = rslt$adj.P.Val
  rank_vector = abs(rslt$t)
  names(rank_vector) = rownames(rslt)
  degs = rownames(rslt[rslt$padj < qval, ])
  return(list(table = rslt, rank_vector = rank_vector, degs = degs))
}

#' @title Helper Function For DEG Analysis of microArray Data (for all cases, including single and multiple comparisons)
#' @param eset Processed gene expression data in the ExpressionSet format
#'              batch & experimental conditions are provided in pData. 
#' @param compList a vector of multiple comparisons in the format of contrasts in limma (e.g. c("A-B", "A-C", "B-C"))
#' @param qval q-value cutoff for DEG analysis (default: 0.05)
#' @return DEresult: a list of DEG results, including those for each single comparison and those for the overall comparison.
#'         Each DEG result is in the format of A list containing:
#'         table: table of DEG results.
#'         rank_vector: a vector of t-statistics for every gene.
#'         degs: a vector of gene symbols for DEGs.
#' @export
#' @import edgeR
#' @import limma
MicroDegs = function(eset, compList, qval = 0.05) {
  
  data = exprs(eset)
  phenodata = pData(eset)
  celltype = as.character(phenodata$celltype)
  
  if (missing(compList) || (length(compList) == 1)) {
    
    return(DEG_Analysis_Micro(eset = eset, qval = qval))
    
  } else {
    
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    
    for (comp in compList) {
      
      cs = strsplit(comp, split = "-")[[1]]
      ids = which(celltype %in% cs)
      tmpdata = data[, ids]
      
      if (is.null(phenodata$batch)){
        tmppdata = new("AnnotatedDataFrame", data = data.frame(celltype = celltype[ids]))
      }else{
        tmppdata = new("AnnotatedDataFrame", data = data.frame(celltype = celltype[ids], batch = batch[ids]))
      }
      
      rownames(tmppdata) = colnames(tmpdata)
      tmpeset = ExpressionSet(assayData = tmpdata, phenoData = tmppdata)
      tmpDErslt = DEG_Analysis_Micro(tmpeset, qval)
      DEresult[[comp]] = tmpDErslt
      
    }
    
    celltype = factor(phenodata$celltype)
    levs = levels(celltype)
    
    if (is.null(phenodata$batch)){
      design = model.matrix(~ 0 + celltype)
      colnames(design) = levs
      contr.matrix = makeContrasts(contrasts = compList, levels = levs)
    }else{
      batch = phenodata$batch
      design = model.matrix(~ 0 + batch + celltype)
      colnames(design) = c("batch", levs)
      contr.matrix = makeContrasts(contrasts = compList, levels = levs)
    }
    vfit = lmFit(data, design)
    vfit = contrasts.fit(vfit, contrasts=contr.matrix)
    efit = eBayes(vfit)
    DEtable = topTable(efit, coef = NULL, number = Inf, sort.by = "F")
    rank_vector = DEtable$F; names(rank_vector) = rownames(DEtable)
    DEtable$padj = DEtable$adj.P.Val
    degs = rownames(DEtable[DEtable$padj < qval, ])
    DEresult$Overall = list(table = DEtable, rank_vector = rank_vector, degs = degs)
    return(DEresult)
    
  }
}

#' @title Helper Function For DEG Analysis of RNA-seq Data using limma + Voom
#' @param counts Processed gene expression count data
#' @param phenodata pData that provides batch & experimental conditions
#' @param complist a vector of multiple comparisons in the format of contrasts in limma (e.g. c("A-B", "A-C", "B-C"))
#' @param lfc (optional) log fold change constraints for DEGs
#' @param qval q-value cutoff for DEG analysis (default: 0.05)
#' @return DEresult: a list of DEG results, including those for each single comparison and those for the overall comparison.
#'         Each DEG result is in the format of A list containing:
#'         table: table of DEG results.
#'         rank_vector: a vector of t-statistics for every gene.
#'         degs: a vector of gene symbols for DEGs.
#'         e: expression data (CPM).
#'         e_batch: batch corrected expression.
#' @export
#' @import edgeR
#' @import limma
RNAseqDegs_limma = function(counts, phenodata, complist, lfc, qval = 0.05) {
#  require(edgeR)
#  require(limma)
  celltype = phenodata$celltype
  levs = levels(celltype)
  dge = DGEList(counts=counts,group=celltype, genes = rownames(counts))
  dge <- calcNormFactors(dge)
  
  if (missing(complist)) {
    
    if(length(levs) > 2) print("Warning ... more than two types!")
    compList = paste(levs, collapse = "-")
    
  } else {
    
    compList <- complist
    
  }
  
  if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
    if_batch = FALSE
    design = model.matrix(~ 0 + celltype)
    colnames(design) = levs
  }else{
    if_batch = TRUE
    batch = phenodata$batch
    levs_batch = levels(batch)
    design = model.matrix(~ 0 + celltype + batch)
    colnames(design) = c(levs, levs_batch[2:length(levs_batch)])
  }
  
  contr.matrix = makeContrasts(contrasts = compList, levels = design)
  v <- voom(dge, design, plot = FALSE)
  vfit = lmFit(v, design)
  vfit = contrasts.fit(vfit, contrasts=contr.matrix)
  efit = eBayes(vfit)
  results = decideTests(efit, p.value = qval)
  
  #log fold change restriction
  if (!missing(lfc)) {
    
    efit <- treat(vfit, lfc = lfc)
    results <- decideTests(efit, p.value = qval)
  }
  print(summary(results))
  
  #remove batch effects
  if(if_batch) {
    design1 = model.matrix(~ 0 + celltype)
    data.rm.batch <- removeBatchEffect(v$E, batch, design=design1)  # limma remove batch
  } else{
    data.rm.batch = NULL
  }
  
  if (length(compList) == 1) {
    
    rslt = topTable(efit, coef=compList, number = Inf, sort.by = "P")
    rslt$padj = rslt$adj.P.Val
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< qval, ])
    
    return(list(table = rslt, rank_vector = rank_vector, degs = degs, 
                e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch)))
    
  } else {
    
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    for (comp in compList){
      rslt = topTable(efit, coef = comp, number = Inf, sort.by = "P")
      rslt$padj = rslt$adj.P.Val
      rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
      degs = rownames(rslt[rslt$padj< qval, ])
      DEresult[[comp]] = list(table = rslt, rank_vector = rank_vector, degs = degs)
    }
    
    if (missing(lfc)) {
      
      DEtable = topTable(efit, coef = NULL, number = Inf, sort.by = "F")
      rank_vector = DEtable$F; names(rank_vector) = rownames(DEtable)
      
    } else {
      
      DEtable = topTable(efit, coef = NULL, number = Inf, sort.by = "P")
      rank_vector = abs(DEtable$t); names(rank_vector) = rownames(DEtable)
    }
    
    DEtable$padj = DEtable$adj.P.Val
    degs = rownames(DEtable[DEtable$padj < qval, ])
    DEresult[["Overall"]] = list(table = DEtable, rank_vector = rank_vector, degs = degs,
                                 e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch))
    
    return(DEresult)
    
  }
}

#' @title Helper Function For DEG Analysis of RNA-seq Data using DESeq
#' @param counts Processed gene expression count data
#' @param phenodata pData that provides batch & experimental conditions
#' @param complist a vector of multiple comparisons in the format of contrasts in limma (e.g. c("A-B", "A-C", "B-C"))
#' @param qval q-value cutoff for DEG analysis (default: 0.05)
#' @return DEresult: a list of DEG results, including those for each single comparison and those for the overall comparison.
#'         Each DEG result is in the format of A list containing:
#'         table: table of DEG results.
#'         rank_vector: a vector of t-statistics for every gene.
#'         degs: a vector of gene symbols for DEGs.
#'         e: expression data (CPM).
#' @export
#' @import DESeq2
RNAseqDegs_DESeq = function(counts, phenodata, complist, qval = 0.05) {
#  require(DESeq2)
  celltype = phenodata$celltype
  
  if (missing(complist) | length(complist) == 1) {
    
    if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
      dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
    }else{
      dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ batch + celltype)
    }
    
    #normalized expression
    dds <- estimateSizeFactors(dds)
    e <- counts(dds, normalized = TRUE)
    
    dds = DESeq(dds)
    DEtable = results(dds, alpha = qval)
    DEtable = data.frame(DEtable)
    DEtable = DEtable[complete.cases(DEtable),]
#    DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
    degs = rownames(DEtable[DEtable$padj < qval, ])
    return(list(table = DEtable, rank_vector = rank_vector, degs = degs, e = e))
    
  } else {
    
    compList <- complist
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    
    for (comp in compList) {
      cs = strsplit(comp, split = "-")[[1]]
      ids = which(celltype %in% cs)
      tmpcounts = counts[, ids]
      tmpph = as(phenodata[ids, ], "data.frame")
      
      if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
        dds = DESeqDataSetFromMatrix(countData = data.frame(tmpcounts), colData = tmpph, design = ~ celltype)
      }else{
        dds = DESeqDataSetFromMatrix(countData = data.frame(tmpcounts), colData = tmpph, design = ~ batch + celltype)
      }
      
      dds = DESeq(dds)
      DEtable = results(dds, alpha = qval)
      DEtable = data.frame(DEtable)
      DEtable = DEtable[complete.cases(DEtable),]
      #DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
      DEtable = DEtable[order(DEtable$pvalue), ]
      rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
      degs = rownames(DEtable[DEtable$padj < qval, ])
      tmpDErslt <- list(table = DEtable, rank_vector = rank_vector, degs = degs)
      DEresult[[comp]] = tmpDErslt
    }
    
    #Multi Comparison
    phenodata = as(phenodata, "data.frame")
    
    if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
      dds = DESeqDataSetFromMatrix(countData = data.frame(counts), colData = phenodata, design = ~ celltype)
    }else{
      dds = DESeqDataSetFromMatrix(countData = data.frame(counts), colData = phenodata, design = ~ batch + celltype)
    }
    
    #normalized expression
    dds <- estimateSizeFactors(dds)
    e <- counts(dds, normalized = TRUE)
    
    #DESeq
    dds = DESeq(dds, test = "LRT", reduced = ~ 1)
    DEtable = results(dds, alpha = qval)
    DEtable = data.frame(DEtable)
    DEtable = DEtable[complete.cases(DEtable),]
#    DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
    degs = rownames(DEtable[DEtable$padj < qval, ])
    
    DEresult$Overall = list(table = DEtable, rank_vector = rank_vector, degs = degs, e = e)
    return(DEresult)
    
  }
  
}
