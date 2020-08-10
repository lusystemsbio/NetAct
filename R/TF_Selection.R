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



GSEA_new = function(GSDB, DElist, minSize=5, nperm = 1000, qval = T){
    stats_vector = as.numeric(DElist)
    gene_list = names(DElist); N = length(gene_list)
    GSDB = lapply(GSDB, intersect, gene_list)
    lens = sapply(GSDB, length)
    ids = which(lens >= minSize)
    GSDB = GSDB[ids]
    gsea_rslt = do.call(rbind, lapply(GSDB, GSEA, gene_list = gene_list, stats_vector= stats_vector))
    # Index the GS
    GSDB = (lapply(GSDB, function(x) fmatch(x, gene_list, nomatch = 0)))
    message("start calculating ......")
    sim_rslt = GSEA_permute(GSDB = GSDB, stats_vector = stats_vector, nperm = nperm)
    message("finish calculating done !!!!!!")
    sim_rslt = cbind(gsea_rslt, sim_rslt)
    pvals = apply(sim_rslt, 1, function(x) sum(x[1] < x)/ nperm)
    z = apply(sim_rslt, 1, function(x) (x[1]-mean(x[2:nperm+1]))/sd(x[2:nperm+1]) )
    rslt = data.frame(tf = rownames(sim_rslt), es = gsea_rslt,
                      lens = lens[ids], pvals = pvals, z = z)
    rslt = rslt[order(rslt$pvals, decreasing = F),]
    if(qval){
        qvals = qvalue(rslt$pvals, lambda=0)$qvalues
        rslt$qvals = qvals
    }
    return(rslt)
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


#####################################################################
#########################  Summary Function #########################
#####################################################################

## GSEA Step
TF_GSEA = function(GSDB, DErslt, minSize=5, nperm = 5000, specific = NULL, qval = T){
    DElist = DErslt$rank_vector
    if (! is.null(specific)){
        stopifnot(specific %in% c("fast", "old"))
        if (specific == "fast"){
            gseaRes = fgsea(GSDB, DElist, nperm=nperm, maxSize=1000, minSize = minSize)
            gseaRes = gseaRes[order(gseaRes$padj),]
            gseaRes = data.frame(gseaRes[,c(1:3, 7)])
            colnames(gseaRes)[1] = "tf"
            if(qval){
                gseaRes$qval = qvalue(gseaRes$pval)$qvalues
            }
            return(gseaRes)
        }else{
            l = permut_glist(GSDB, DElist, minSize = minSize, n_permutation = nperm)
            return(l)
        }
    }else{
        gseaRes = GSEA_new(GSDB, DElist, minSize = minSize, nperm = nperm, qval = qval)
        return(gseaRes)
    }
}


## TF selection step:
TF_Selection_old = function(GSEArslt, DErslt, qval = 0.05, GSDB) {
    table = DErslt$table
    degs = rownames(table[table$padj< 0.05, ])
    tfs = GSEArslt$tf[which(GSEArslt$qval < qval)]
    c = cal_coverage(tfs, GSDB, degs)
    cat(length(tfs)," TFs found!\n")
    cat(TF_list_IPA(tfs),"\n")
    return(list(tfs=tfs, c=c, degs = degs))
}


cal_coverage = function(tfs, GSDB, degs) {
    coverage = numeric(); All_targets = character()
    ndegs = length(degs)
    for (tf in tfs) {
        All_targets = unique(c(All_targets, GSDB[[tf]]))
        tmp_coverage = sum(All_targets %in% degs)/ ndegs
        coverage=c(coverage, tmp_coverage)
    }
    names(coverage) = tfs
    plot(coverage, col="red", type="o",pch=18)
    return(coverage)
}

TF_list_IPA = function(tfs) {
    return(paste(sort(tfs), collapse = ","))
}

## TF selection step:
TF_Selection = function(GSEArslt,qval = 0.05, compList = NULL, ntop = NULL) {
  if(is.null(compList) | length(compList) == 1) {
    tfs = TF_Selection_single(GSEArslt, qval = qval, ntop = ntop)
  }else{
    tfs = character()
    i = 0
    for(comp in compList){
      i = i + 1
      tfs_single = TF_Selection_single(GSEArslt[[i]], qval = qval, ntop = ntop)
      tfs = c(tfs, tfs_single)
    }
    tfs = sort(unique(tfs))
  }
  return(tfs)
}

TF_Selection_single = function(GSEArslt, qval = 0.05, ntop = NULL) {
  data_rslt = GSEArslt[order(GSEArslt[,7], GSEArslt[,6], decreasing = c(FALSE, TRUE)),]

  tfs = as.character(data_rslt$tf[data_rslt$qvals < qval])
  if(!is.null(ntop)){
    if(length(tfs) > ntop){
      tfs = tfs[1:ntop]
    }
  }
  tfs = sort(tfs)
  return(tfs)
}


