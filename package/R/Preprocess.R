##################################################################
################### DE gene analysis functions ###################
##################################################################

MicroDegs = function(eset){
    data = exprs(eset)
    #the first column of the phenotype data is the celltype
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
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< 0.05, ])
    return(list(table = rslt, rank_vector = rank_vector, degs = degs))
}

## Gene ranking list from RNAseq data
RNAseqDegs = function(counts, phenodata){
    require(DESeq2)
    if (is.null(phenodata$batch) | length(unique(phenodata$batch)) == 1){
        dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
    }else{
        dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ batch + celltype)
    }
    dds = DESeq(dds)
    DEtable = results(dds)
    DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
    return(list(table = DEtable, rank_vector = rank_vector))
}

## Multiple Group Comparison to get the DE result 
MultiMicroDegs = function(eset, compList){
    # compList includes the comparison names e.g. compList = c("Early-Middle", "Middle-Late", "Early-Late")
    data = exprs(eset)
    phenodata = pData(eset)
    celltype = phenodata$celltype
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
    degs = rownames(DEtable[DEtable$padj< 0.05, ])
    return(list(table = DEtable, rank_vector = rank_vector, degs = degs))
}

MultiRNAseqDegs = function(counts, phenodata){
    require(DESeq2)
    if (is.null(phenodata$batch)){
        dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
    }else{
        dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ batch + celltype)
    }
    dds = DESeq(dds, test = "LRT", reduced = ~ 1)
    DEtable = results(dds)
    DEtable =  data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
    return(list(table = DEtable, rank_vector = rank_vector))
}


# preprocess data
MicroProcess = function(eset, compList) {
    data = exprs(eset)
    #the first column of the phenotype data is the celltype
    phenodata = pData(eset)
    celltype = as.character(phenodata$celltype)
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    for (comp in compList){
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
        tmpDErslt = MicroDegs(tmpeset)
        DEresult[[comp]] = tmpDErslt
    }
    DEresult$Overall = MultiMicroDegs(eset, compList)
    return(DEresult)
}

# preprocess data
RNAseqProcess = function(counts, phenodata, compList) {
    require(DESeq2)
    celltype = phenodata$celltype
    DEresult = vector(mode = "list", length = length(compList) + 1)
    names(DEresult) = c(compList, "Overall")
    for (comp in compList){
      cs = strsplit(comp, split = "-")[[1]]
      ids = which(celltype %in% cs)
      tmpcounts = counts[, ids]
      tmpph = phenodata[ids, ]
      tmpDErslt = RNAseqDegs(tmpcounts, tmpph)
      DEresult[[comp]] = tmpDErslt
    }
    DEresult$Overall = MultiRNAseqDegs(counts, phenodata)
    return(DEresult)
}

