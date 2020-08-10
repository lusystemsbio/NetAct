##################################################################
################### DE gene analysis functions ###################
##################################################################

MicroDegs = function(eset, qval = 0.05){
    require(limma)
    data = exprs(eset)
    #the first column of the phenotype data is the celltype
    phenodata = pData(eset)
    celltype = phenodata$celltype
    levs = levels(celltype)
    if(length(levs) > 2) print("Warning ... more than two types!")
    compList = paste(levs, collapse = "-")
    
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
    fit = lmFit(data, design)
    fit = contrasts.fit(fit, contrasts=contr.matrix)
    efit = eBayes(fit)
    results = decideTests(efit, p.value = qval)
    print(summary(results))
    
    rslt = topTable(efit, coef = compList, number = Inf, sort.by = "P")
    rslt$padj = rslt$adj.P.Val
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< qval, ])
    
    if(if_batch) {
      design1 = model.matrix(~ 0 + celltype)
      data.rm.batch <- removeBatchEffect(data, batch, design=design1)  # limma remove batch
    } else{
      data.rm.batch = NULL
    }

    return(list(table = rslt, rank_vector = rank_vector, degs = degs, 
                e_batch = as.data.frame(data.rm.batch)))
}

RNAseqDegs_limma = function(counts, phenodata, qval = 0.05){
    require(edgeR)
    require(limma)
    celltype = phenodata$celltype
    levs = levels(celltype)
    if(length(levs) > 2) print("Warning ... more than two types!")
    compList = paste(levs, collapse = "-")

    dge = DGEList(counts=counts,group=celltype, genes = rownames(counts))
    dge <- calcNormFactors(dge)
    
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
    v <- voom(dge, design, plot=TRUE)
    vfit = lmFit(v, design)
    vfit = contrasts.fit(vfit, contrasts=contr.matrix)
    efit = eBayes(vfit)
    results = decideTests(efit, p.value = qval)
    print(summary(results))
    rslt = topTable(efit, coef=compList, number = Inf, sort.by = "P")
    rslt$padj = rslt$adj.P.Val
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< qval, ])
    
    if (if_batch) {
      design1 = model.matrix(~ 0 + celltype)
      data.rm.batch <- removeBatchEffect(v$E, batch, design=design1)  # limma remove batch
    }else{
      data.rm.batch = NULL
    }
    
    return(list(table = rslt, rank_vector = rank_vector, degs = degs, 
                e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch)))
}

## Multiple Group Comparison to get the DE result 
MultiMicroDegs = function(eset, compList, qval = 0.05){
    # compList includes the comparison names e.g. compList = c("Early-Middle", "Middle-Late", "Early-Late")
    data = exprs(eset)
    phenodata = pData(eset)
    celltype = phenodata$celltype
    levs = levels(celltype)
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
    fit = lmFit(data, design)
    fit = contrasts.fit(fit, contrasts=contr.matrix)
    efit = eBayes(fit)
    results = decideTests(efit, p.value = qval)
    print(summary(results))
    
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    for (comp in compList){
      rslt = topTable(efit, coef = comp, number = Inf, sort.by = "P")
      rslt$padj = rslt$adj.P.Val
      rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
      degs = rownames(rslt[rslt$padj< qval, ])
      DEresult[[comp]] = list(table = rslt, rank_vector = rank_vector, degs = degs)
    }
    
    if (if_batch) {
      design1 = model.matrix(~ 0 + celltype)
      data.rm.batch <- removeBatchEffect(data, batch, design=design1)  # limma remove batch
    }else{
      data.rm.batch = NULL
    }
    
    DEtable = topTable(efit, coef = NULL, number = Inf, sort.by = "F")
    rank_vector = DEtable$F; names(rank_vector) = rownames(DEtable)
    DEtable$padj = DEtable$adj.P.Val
    degs = rownames(DEtable[DEtable$padj< qval, ])
    DEresult[["Overall"]] = list(table = DEtable, rank_vector = rank_vector, degs = degs, 
                                 e_batch = as.data.frame(data.rm.batch))
    
    return(DEresult)
}

MultiRNAseqDegs_limma = function(counts, phenodata, compList, qval = 0.05){
  require(edgeR)
  require(limma)
  celltype = phenodata$celltype
  levs = levels(celltype)
  dge = DGEList(counts=counts,group=celltype, genes = rownames(counts))
  dge <- calcNormFactors(dge)
  
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
  v <- voom(dge, design, plot=TRUE)
  vfit = lmFit(v, design)
  vfit = contrasts.fit(vfit, contrasts=contr.matrix)
  efit = eBayes(vfit)
  results = decideTests(efit, p.value = qval)
  print(summary(results))
  
  DEresult = vector(mode = "list", length = length(compList) + 1)
  names(DEresult) = c(compList, "Overall")
  for (comp in compList){
    rslt = topTable(efit, coef = comp, number = Inf, sort.by = "P")
    rslt$padj = rslt$adj.P.Val
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< qval, ])
    DEresult[[comp]] = list(table = rslt, rank_vector = rank_vector, degs = degs)
  }
  
  if(if_batch) {
    design1 = model.matrix(~ 0 + celltype)
    data.rm.batch <- removeBatchEffect(v$E, batch, design=design1)  # limma remove batch
  } else{
    data.rm.batch = NULL
  }
  
  DEtable = topTable(efit, coef = NULL, number = Inf, sort.by = "F")
  rank_vector = DEtable$F; names(rank_vector) = rownames(DEtable)
  DEtable$padj = DEtable$adj.P.Val
  degs = rownames(DEtable[DEtable$padj< qval, ])
  DEresult[["Overall"]] = list(table = DEtable, rank_vector = rank_vector, degs = degs,
                               e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch))
  
  return(DEresult)
}

## Gene ranking list from RNAseq data
RNAseqDegs_DEseq = function(counts, phenodata, qval = 0.05){
  require(DESeq2)
  if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
    dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
  }else{
    dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ batch + celltype)
  }
  dds = DESeq(dds)
  DEtable = results(dds, alpha = qval)
  DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
  DEtable = DEtable[order(DEtable$pvalue), ]
  rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
  return(list(table = DEtable, rank_vector = rank_vector))
}

MultiRNAseqDegs_DEseq = function(counts, phenodata, qval = 0.05){
    require(DESeq2)
    if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
        dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
    }else{
        dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ batch + celltype)
    }
    dds = DESeq(dds, test = "LRT", reduced = ~ 1)
    DEtable = results(dds, alpha = qval)
    DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
    return(list(table = DEtable, rank_vector = rank_vector))
}

# preprocess data
RNAseqProcess_DEseq = function(counts, phenodata, compList) {
    require(DESeq2)
    celltype = phenodata$celltype
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    for (comp in compList){
      cs = strsplit(comp, split = "-")[[1]]
      ids = which(celltype %in% cs)
      tmpcounts = counts[, ids]
      tmpph = phenodata[ids, ]
      tmpDErslt = RNAseqDegs_DEseq(tmpcounts, tmpph)
      DEresult[[comp]] = tmpDErslt
    }
    DEresult$Overall = MultiRNAseqDegs_DEseq(counts, phenodata)
    return(DEresult)
}
