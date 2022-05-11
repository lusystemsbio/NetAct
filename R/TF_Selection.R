##########################################################
############# TF GSEA and Selection Function #############
##########################################################
## GSEA implementation
GSEA = function(gene_list, gene_set, stats_vector){
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
GSEA_permutone= function(sim_all, gene_set, stats_vector, N){
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

####################################################################
###################### Old permutation method ######################
####################################################################
## Permutate once
permut_one = function(sim_all, gs, stats_vector){
  tmp_sim_gseas = unlist(lapply(sim_all, function(x) GSEA(x, gs, stats_vector)))
  return(tmp_sim_gseas)
}

## Normal mixture model p value
nm_pval=function(permutated_vector){
  ts=permutated_vector[1]
  tmp_m=Mclust(permutated_vector[2:length(permutated_vector)], parameter = T, modelNames = "V")
  pval=0
  for (i in 1:length(tmp_m$parameters$mean)){
    tmp_prob=as.numeric(tmp_m$parameters$pro[i])
    tmp_mean=as.numeric(tmp_m$parameters$mean[i])
    tmp_sd=sqrt(as.numeric(tmp_m$parameters$variance$sigmasq[i]))
    pval=pval+ tmp_prob* (1-pnorm(ts, tmp_mean, tmp_sd))
  }
  return(pval)
}

## Permutation of the list
permut_glist = function(gs_db, gene_list, minSize=5, n_permutation = 200){
  genes = names(gene_list)
  n_lens = sapply(gs_db, function(x) length(intersect(x, genes)))
  ids = which(n_lens >= minSize)
  gs_db = gs_db[ids]
  sim_glists = vector(mode = "list", n_permutation)
  sim_glists = lapply(sim_glists, function(x) x = sample(names(gene_list)))
  sim_all = append(list(names(gene_list)), sim_glists)
  tmp_stats_vector = as.numeric(gene_list)
  
  rlst_vector = numeric()
  for (i in 1:length(gs_db)){
    tmp_sim_gseas = permut_one(sim_all, gs_db[[i]], tmp_stats_vector)
    tmp_gsea = tmp_sim_gseas[1]
    if (tmp_gsea == -1){
      rlst_vector = c(rlst_vector, -1, 1, 1)
    }
    else{
      cat(names(gs_db)[i], "\n")
      pval1 = nm_pval(tmp_sim_gseas)
      pval2 = sum(tmp_sim_gseas > tmp_sim_gseas[1])/ n_permutation
      rlst_vector = c(rlst_vector,tmp_gsea, pval1, pval2)
    }
  }
  rslt_mat = matrix(rlst_vector, ncol = 3, byrow = T)
  rslt_mat = cbind(rslt_mat, qvalue(rslt_mat[,2])$qvalues,  qvalue(rslt_mat[,3])$qvalues, n_lens[ids])
  colnames(rslt_mat) =  c("gsea", "pval1", "pval2", "qval1", "qval2", "size")
  rslt_mat = rslt_mat[order(rslt_mat[,2], decreasing = F),]
  return(data.frame(rslt_mat))
}

###
GSEA_new = function(GSDB, DElist, minSize=5, nperm = 1000) {
  stats_vector = as.numeric(DElist)
  gene_list = names(DElist); N = length(gene_list)
  GSDB = lapply(GSDB, intersect, gene_list)
  lens = sapply(GSDB, length)
  ids = which(lens >= minSize)
  GSDB = GSDB[ids]
  gsea_rslt = do.call(rbind, lapply(GSDB, NetAct::GSEA, gene_list = gene_list, stats_vector= stats_vector))
  # Index the GS
  GSDB = (lapply(GSDB, function(x) fmatch(x, gene_list, nomatch = 0)))
  message("start calculating ......")
  sim_rslt = NetAct::GSEA_permute(GSDB = GSDB, stats_vector = stats_vector, nperm = nperm)
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


## GSEA Step
TF_GSEA = function(GSDB, DErslt, minSize=5, nperm = 5000, specific = NULL){
  DElist = DErslt$rank_vector
  if (! is.null(specific)){
    stopifnot(specific %in% c("fast", "old"))
    if (specific == "fast"){
      gseaRes = fgsea(GSDB, DElist, nperm=nperm, maxSize=1000, minSize = minSize)
      gseaRes = gseaRes[order(gseaRes$padj),]
      gseaRes = data.frame(gseaRes[,c(1:3, 7)])
      colnames(gseaRes)[1] = "tf"
      gseaRes$qval = qvalue(gseaRes$pval)$qvalues
      
      return(gseaRes)
    }else{
      l = permut_glist(GSDB, DElist, minSize = minSize, n_permutation = nperm)
      return(l)
    }
  }else{
    gseaRes = GSEA_new(GSDB, DElist, minSize = minSize, nperm = nperm)
    return(gseaRes)
  }
}



####TF_Selection####
TF_Selection = function(GSDB, DErslt, minSize=5, nperm = 5000, specific = NULL, qval = 0.05, compList = NULL, ntop = NULL, nameFile = NULL) {
  
  if(is.null(compList) | length(compList) == 1) {
    
    gsearslt <- TF_GSEA(GSDB = GSDB, DErslt = DErslt, minSize = minSize, nperm = nperm, specific = specific)
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
      
      gsea_results_list[[i]] <- TF_GSEA(GSDB = GSDB, DErslt = DErslt[[i]], minSize = minSize, nperm = nperm, specific = specific)
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

####Reselect_TFs####
Reselect_TFs = function(GSEArslt, qval = 0.05, combine_TFs = TRUE, ntop = NULL) {
  
  if (length(GSEArslt) == 1) {
    
    tfs <- as.character(GSEArslt$tf[GSEArslt$qvals < qval[1]])
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
