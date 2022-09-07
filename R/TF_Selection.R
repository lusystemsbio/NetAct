##########################################################
############# TF GSEA and Selection Function #############
##########################################################

#' @title Compute the enrichment score (ES) from Gene Set Enrichment Analysis
#' @description In this specific implementation, we modified two aspects:
#'        (1) the stats_vector takes the absolute value of the t-statistics; 
#'        This ensures that the GSEA works for the case where part of genes are up-regulated and others are down-regulated.
#'        (2) permutation of the gene symbols of the ranking vector, instead of gene expression values.
#' @param gene_list a vector of genes in the expression data
#' @param gene_set a vector of genes in the gene set
#' @param stats_vector a vector of DEG statistics for every gene in gene_list (rank_vector in the DEG results)
#'                     Absolute values of the t-statistics are required for the desired performance.
#' @return ES: enrichment score
#' @export
GSEA_score = function(gene_list, gene_set, stats_vector){
  # Input genelist must be ordered.
  tag_indicator = sign(fmatch(gene_list, gene_set, nomatch = 0))
  no_tag_indicator = 1- tag_indicator
  N = length(gene_list); Nh = length(gene_set); Nm =  N - Nh
  
  sum_rank_tag = sum(stats_vector[tag_indicator ==1])
  if (sum_rank_tag == 0){
    ES = -1
  }else{
    norm_tag = 1.0/ sum_rank_tag; norm_no_tag = 1/Nm
    RES = cumsum(tag_indicator*stats_vector*norm_tag - no_tag_indicator*norm_no_tag)
    max_ES = max(RES)
    min_ES = min(RES)
    if (max_ES > - min_ES) {
      ES = max_ES
    } else {
      ES = min_ES
    }
  }
  return(ES)
}

####################################################################
###################### New permutation method ######################
####################################################################

#' @title Compute ES scores for Gene Set Enrichment Analysis (GSEA) with a new permutation method (using a revised algorithm)
#' @description To improve computational efficiency, we devised a new permutation approach by swapping stats_vector. 
#'              Here, the gene symbols/names are permutated without changing the ranking vector (stats_vector).
#'              This function becomes unused in NetAct, as a much faster c++ implementation (GSEA_permute) is provided.
#' @param sim_all a vector of genes in the expression data
#' @param gene_set a vector of genes in the gene set
#' @param stats_vector a vector of DEG statistics for every gene in gene_list (rank_vector in the DEG results);
#'                     Absolute values of the t-statistics are required for the desired performance.
#' @param N total number of genes (size of sim_all)
#' @return ES: enrichment score
GSEA_permut_R_revised = function(sim_all, gene_set, stats_vector, N){
  Nh = length(gene_set); Nm =  N - Nh; norm_no_tag = 1/Nm
  ES_cal_mat = lapply(sim_all, function(gene_list){
    random_tag_indicator = sign(fmatch(gene_list, gene_set, nomatch = 0))
    random_no_tag_indicator = 1 - random_tag_indicator
    random_sum_rank_tag = sum(stats_vector[random_tag_indicator ==1])
    random_norm_tag = 1.0/ random_sum_rank_tag
    random_RES_vec = random_tag_indicator*stats_vector*random_norm_tag -
      random_no_tag_indicator*norm_no_tag
    return(random_RES_vec)
  })
  ES_cal_mat = do.call(cbind, ES_cal_mat)
  RES_mat = apply_cumsum_col(ES_cal_mat)
  ES_vec = apply(RES_mat, 2, function(RES){
    max_ES = max(RES)
    min_ES = min(RES)
    ifelse(max_ES > -min_ES, max_ES, min_ES)
  })
  ES_vec = unlist(ES_vec)
  return(ES_vec)
}

#' @title  Compute ES scores for Gene Set Enrichment Analysis (GSEA) with a new permutation method (using the original GSEA algorithm)
#' @description The function uses the original GSEA enrichment score calculation but using the new permutation method.
#'              Here, the gene symbols/names are permutated without changing the ranking vector (stats_vector).
#' @param sim_all a matrix of permutated gene lists
#' @param gs a vector of genes in the gene set
#' @param stats_vector a vector of DEG statistics for every gene in gene_list (rank_vector in the DEG results)
#'                     Absolute values of the t-statistics are required for the desired performance.
#' @return tmp_sim_sgeas: a vector of ES values for all permutated gene lists
GSEA_permut_R = function(sim_all, gs, stats_vector){
  tmp_sim_gseas = unlist(lapply(sim_all, function(x) GSEA_score(x, gs, stats_vector)))
  return(tmp_sim_gseas)
}

#' @title Gene Set Enrichment Analysis (GSEA) with a new permutation method -- implementation in R
#' @description To improve computational efficiency, we devised a new permutation approach by swapping stats_vector. 
#'              Here, the gene symbols/names are permutated without changing the ranking vector (stats_vector).
#' @param GSDB gene set database (a list of gene sets, each of which is comprised of a vector genes)
#' @param DElist a vector of DEG statistics for every gene in gene_list (rank_vector in the DEG results)
#' @param minSize the minimum number of overlapping genes required for each gene set (a gene set filtering parameter, default: 5)
#' @param nperm the number of gene list permutations (default: 1000)
#' @return data.frame(rslt_mat): a table of GSEA results:
#'         tf: TF (gene set name).
#'         es: ES score.
#'         lens: number of overlapping genes in each gene set.
#'         pvals: p-value by counting.
#'         z: z-score.
#'         qvals: q-value from pvals.
#' @import qvalue 
#' @export
GSEA_proc_R = function(GSDB, DElist, minSize=5, nperm = 1000){
  
  genes = names(DElist)
  n_lens = sapply(GSDB, function(x) length(intersect(x, genes)))
  ids = which(n_lens >= minSize)
  GSDB = GSDB[ids]
  sim_glists = vector(mode = "list", nperm)
  sim_glists = lapply(sim_glists, function(x) x = sample(names(DElist)))
  sim_all = append(list(names(DElist)), sim_glists)
  tmp_stats_vector = as.numeric(DElist)
  
  n_DB = length(GSDB)
  sim_rslt = matrix(0, nrow = n_DB , ncol = nperm + 1)
  for (i in 1:length(GSDB)){
    sim_rslt[i,] = GSEA_permut_R(sim_all, GSDB[[i]], tmp_stats_vector)
  }
  rownames(sim_rslt) = names(GSDB)
  
  pvals = apply(sim_rslt, 1, function(x) sum(x[1] < x)/ nperm)
  z = apply(sim_rslt, 1, function(x) (x[1]-mean(x[2:nperm+1]))/sd(x[2:nperm+1]) )
  rslt = data.frame(tf = rownames(sim_rslt), es = sim_rslt[,1],
                    lens = n_lens[ids], pvals = pvals, z = z)
  rslt = rslt[order(rslt$pvals, decreasing = F),]
  
  qvals = qvalue(rslt$pvals, lambda=0)$qvalues
  rslt$qvals = qvals
  
  return(rslt)
}

#' @title Gene Set Enrichment Analysis (GSEA) with a new permutation method -- implementation in R/c++ 
#' @description To improve computational efficiency, we devised a new permutation approach by swapping stats_vector. 
#'              Here, the gene symbols/names are permutated without changing the ranking vector (stats_vector).
#'              A much faster c++ implementation (GSEA_permute) is used.
#' @param GSDB gene set database (a list of gene sets, each of which is comprised of a vector genes)
#' @param DElist a vector of DEG statistics for every gene in gene_list (rank_vector in the DEG results)
#' @param minSize the minimum number of overlapping genes required for each gene set (a gene set filtering parameter, default: 5)
#' @param nperm the number of gene list permutations (default: 1000)
#' @return data.frame(rslt_mat): a table of GSEA results:
#'         tf: TF (gene set name).
#'         es: ES score.
#'         lens: number of overlapping genes in each gene set.
#'         pvals: p-value by counting.
#'         z: z-score.
#'         qvals: q-value from pvals.
#' @import Rcpp
#' @import qvalue 
#' @import fastmatch
#' @export
GSEA_proc_RC = function(GSDB, DElist, minSize=5, nperm = 1000) {
  stats_vector = as.numeric(DElist)
  gene_list = names(DElist); N = length(gene_list)
  GSDB = lapply(GSDB, intersect, gene_list)
  lens = sapply(GSDB, length)
  ids = which(lens >= minSize)
  GSDB = GSDB[ids]
  gsea_rslt = do.call(rbind, lapply(GSDB, GSEA_score, gene_list = gene_list, stats_vector= stats_vector))
  # Index the GS
  GSDB = (lapply(GSDB, function(x) fmatch(x, gene_list, nomatch = 0)))
  message("start calculating ......")
  sim_rslt = GSEA_permute(GSDB = GSDB, stats_vector = stats_vector, nperm = nperm)
  message("calculations finished !!!!!!")
  sim_rslt = cbind(gsea_rslt, sim_rslt)
  pvals = apply(sim_rslt, 1, function(x) sum(x[1] < x)/ nperm)
  z = apply(sim_rslt, 1, function(x) (x[1]-mean(x[2:nperm+1]))/sd(x[2:nperm+1]) )
  rslt = data.frame(tf = rownames(sim_rslt), es = gsea_rslt,
                    lens = lens[ids], pvals = pvals, z = z)
  rslt = rslt[order(rslt$pvals, decreasing = F),]
  
  qvals = qvalue(rslt$pvals, lambda=0)$qvalues
  rslt$qvals = qvals
  
  return(rslt)
}

#' A unified Gene Set Enrichment Analysis (GSEA) function for three methods
#' @param GSDB gene set database (a list of gene sets, each of which is comprised of a vector genes)
#' @param DErslt DEG results
#' @param minSize the minimum number of overlapping genes required for each gene set (a gene set filtering parameter, default: 5)
#' @param nperm the number of gene list permutations (default: 1000)
#' @param method fast: fgsea; r: R implementation of GSEA with a new permutation method; binary: R/C++ implementation for fast speed
#' @return gseaRes: a table of GSEA results:
#'         tf: TF (gene set name).
#'         es: ES score.
#'         lens: number of overlapping genes in each gene set.
#'         pvals: p-value by counting.
#'         z: z-score.
#'         qvals: q-value from pvals.
#' @import fgsea
#' @import qvalue 
#' @export
TF_GSEA = function(GSDB, DErslt, minSize=5, nperm = 1000, method = "binary"){
  DElist = DErslt$rank_vector
  if(method == "fast"){
    gseaRes = fgsea(GSDB, DElist, nperm=nperm, maxSize=1000, minSize = minSize)
    gseaRes = gseaRes[order(gseaRes$padj),]
    gseaRes = data.frame(gseaRes[,c(1,4,7,2)])
    colnames(gseaRes) = c("tf", "es", "lens", "pvals")
    gseaRes$qvals = qvalue(gseaRes$pval)$qvalues
  }else if(method == "r"){
    gseaRes = GSEA_proc_R(GSDB, DElist, minSize = minSize, nperm = nperm)
  }else{
    gseaRes = GSEA_proc_RC(GSDB, DElist, minSize = minSize, nperm = nperm)
  }
  return(gseaRes)
}

#' Identifying enriched TFs using Gene Set Enrichment Analysis (GSEA) -- a wrapper function with many options
#' @param GSDB gene set database (a list of gene sets, each of which is comprised of a vector genes)
#' @param DErslt DEG results
#' @param minSize the minimum number of overlapping genes required for each gene set (a gene set filtering parameter, default: 5)
#' @param nperm the number of gene list permutations (default: 1000)
#' @param method fast: fgsea; R: r implementation of GSEA with a new permutation method; binary: R/C++ implementation for fast speed
#' @param qval q-value cutoff (default: 0.05)
#' @param compList a vector of comparisons, it needs to be consistent with DErslt from MicroDegs, RNAseqDegs_limma, and RNAseqDegs_DESeq.
#'                  GSEA is applied to each comparison 
#' @param ntop the number of top genes (selection by the top genes) (default: NULL, no selection by the top genes)
#' @param nameFile file name to save the GSEA results (default: NULL, no output to a file). 
#'                  The saved results can be reused later to adjust the TF selection parameters
#' @return a list of results: 
#'         GSEArslt: a dataframe of GSEA results (see TF_GSEA).
#'         tfs: a vector of selected TFs.
#' @export
TF_Selection = function(GSDB, DErslt, minSize=5, nperm = 5000, method = "binary", qval = 0.05, 
                        compList = NULL, ntop = NULL, nameFile = NULL) {
  
  if(is.null(compList) | length(compList) == 1) {
    gsearslt <- TF_GSEA(GSDB = GSDB, DErslt = DErslt, minSize = minSize, nperm = nperm, method = method)
    tfs <- as.character(gsearslt$tf[gsearslt$qvals < qval[1]])
    if(!is.null(ntop)){
      if(length(tfs) > ntop){
        tfs = tfs[1:ntop]
      }
    }
    tfs = sort(tfs)
    output <- list(GSEArslt = list(gsearslt), tfs = tfs)
  }else{
    if ((length(qval) != 1) & (length(qval) != length(compList)))  {
      stop("qvalue length must be equal to 1 OR the number of comparisons")
    }
    if (length(qval) == 1) {
      qval <- rep(qval, length(compList))
    }
    gsea_results_list <- vector(mode = "list", length = length(compList))
    names(gsea_results_list) <- compList
    tfs = character()
    for (i in 1:length(compList)) {
      gsea_results_list[[i]] <- TF_GSEA(GSDB = GSDB, DErslt = DErslt[[i]], minSize = minSize, nperm = nperm, method = method)
      tfs_tmp <- as.character(gsea_results_list[[i]]$tf[gsea_results_list[[i]]$qvals < qval[i]])
      if(!is.null(ntop)){
        if(length(tfs_tmp) > ntop){
          tfs_tmp = tfs_tmp[1:ntop]
        }
      }
      tfs = c(tfs, tfs_tmp)
    }
    tfs = sort(unique(tfs))
    output <- list(GSEArslt = gsea_results_list, tfs = tfs)
  }

  if (!is.null(nameFile)) {
    saveRDS(output, file = paste0(nameFile, ".RDS" ))
  }
  return(output)
}

#' Reselecting TFs using gene set enrichement analysis (GSEA) using an adjusted set of parameters (work together with TF_Selection)
#' @param GSEArslt GSEA results from TF_Selection
#' @param qval q-value cutoff (default: 0.05)
#' @param combine_TFs whether combine selected TFs from multiple comparisons or not (default: TRUE)
#' @param ntop the number of top genes (selection by the top genes) (default: NULL, no selection by the top genes)
#' @return tfs: a vector of selected TFs
#' @export
Reselect_TFs = function(GSEArslt, qval = 0.05, combine_TFs = TRUE, ntop = NULL) {
  if (length(GSEArslt) == 1) {
    tfs <- as.character(GSEArslt[[1]]$tf[which(GSEArslt[[1]]$qvals < qval[1])])
    if(!is.null(ntop)){
      if(length(tfs) > ntop){
        tfs = tfs[1:ntop]
      }
    }
    tfs = sort(tfs)
  } else {
    if ((length(qval) != 1) & (length(qval) != length(GSEArslt)))  {
      stop("qvalue length must be equal to 1 OR the number of comparisons")
    }
    if (length(qval) == 1) {
      qval <- rep(qval, length(GSEArslt))
    }
    tfs = vector(mode = "list", length = length(GSEArslt))
    for (i in 1:length(GSEArslt)) {
      tfs_tmp <- as.character(GSEArslt[[i]]$tf[GSEArslt[[i]]$qvals < qval[i]])
      if(!is.null(ntop)){
        if(length(tfs_tmp) > ntop){
          tfs_tmp = tfs_tmp[1:ntop]
        }
      }
      tfs[[i]] = sort(tfs_tmp)
    }
    if (combine_TFs) {
      tfs = sort(unique(unlist(tfs)))
    } else {
      names(tfs) <- names(GSEArslt)
    }
  }
  return(tfs)
}
