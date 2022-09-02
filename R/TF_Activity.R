########################################################
################ TF Activities Function ################
########################################################

#' Inference of TF activity
#' @param tfs a vector of selected tfs
#' @param GSDB gene set database (a list of gene sets, each of which is comprised of a vector genes)
#' @param eset expression set of gene expression data or gene expression matrix
#' @param DErslt DEG results
#' @param with_weight whether weighting factors (based on DEG p-values) are used to compute TF activity (default: TRUE)
#' @param if_module whether the grouping scheme (activation or inhibition) depends on module detection algorithm (default: FALSE, no need to change)
#' @param ind Hill coefficient parameter used in the weighting factors (default: 1/5, recommend to use 0 < ind < 1/4)
#' @param useCorSign allow the option to use the TF gene expression correlation to flip signs (default: TRUE)
#' @return a list of results: 
#'         all_list: grouping scheme of all TF gene sets.
#'         all_activity: matrix of TF activity.
#' @export
TF_Activity = function (tfs, GSDB, eset, DErslt, with_weight = TRUE, if_module = FALSE,
                          ind = 1/5, useCorSign = TRUE){
    if (is(eset, "ExpressionSet")){
        data = exprs(eset)
    }else{
      data = eset
      }
    table = DErslt$table
    DEweights = data.frame(p = table$padj)
    rownames(DEweights) = rownames(table)
    tfs = tfs[tfs %in% rownames(data)]
    all_activity = all_list = NULL
    
    for (tf in tfs) {
        comms = intersect(GSDB[[tf]], rownames(table))
        g_inds = which(rownames(data) %in% comms)
        tmp_data = data[g_inds, ]
        tmp_data = rem_data(tmp_data)
        cor_mat = cor(as.matrix(t(tmp_data)), method = "spearman")
        tmp_dist <- (cor_mat+1)/2
        diag(tmp_dist) <- 0 
        dimnames(cor_mat) = list(rownames(tmp_data), rownames(tmp_data))
        tmp_names = list(rownames(tmp_data), "sign")
        col_sums = apply(tmp_dist, 2, sum)
        mod_mat = tmp_dist - sapply(col_sums, function(x) x *
                                        col_sums/sum(tmp_dist))
        ev = eigen(mod_mat)
        tmp_sign = data.frame(sign(ev$vectors[, 1]))
        dimnames(tmp_sign) = tmp_names

        ng_mat = cor_mat[tmp_sign$sign == -1, tmp_sign$sign ==
                             -1]
        ps_mat = cor_mat[tmp_sign$sign == 1, tmp_sign$sign ==
                             1]
        test_mats = list(ng_mat = ng_mat, ps_mat = ps_mat)
        gs_list1 = list()
        test_center = numeric()
        for (i in 1:2) {
            test_mat = test_mats[[i]]
            nrow_mat = nrow(test_mat)
            if (is.null(nrow_mat)) {
                next
            }
            else if (nrow_mat <= 2) {
                gs_list1[[i]] = rownames(test_mat)
            }
            else {
                tmp_kmean = kmeans(test_mat, 1)
                centers = tmp_kmean$centers
                test_center = c(test_center, centers)
                test_dim = length(centers)
                distances = apply(test_mat, 1, function(x) cor(as.numeric(x),
                                                               as.numeric(centers), method = "spearman"))
                m = mean(distances)
                names_left = names(distances)[distances >= m]
                gs_list1[[i]] = names_left
            }
        }
        tf_exprs = as.numeric(data[tf, ])
        gs_remain = unlist(gs_list1)
        
        tmp_rslt = cal_activity(gs_remain, tmp_data, tmp_sign,
                                ind, with_weight, DEweights, tf_exprs, useCorSign)
        
        tmp_activity = tmp_rslt$activity
        tmp_sign = tmp_rslt$sign
        gs_list2 = list()
        gs_remain2 = c()
        cors = apply(tmp_data, 1, function(x) cor(x, tf_exprs,
                                                  method = "spearman"))
        pos_cors = cors[cors > 0]
        neg_cors = cors[cors < 0]
        tmp_sign2 = data.frame(sign(cors))
        m_pos = mean(pos_cors)
        m_neg = mean(neg_cors)
        gs_list2[[1]] = names(pos_cors[pos_cors >= m_pos])
        gs_list2[[2]] = names(neg_cors[neg_cors <= m_neg])
        gs_remain = unlist(gs_list2)
        
        tmp_rslt = cal_activity(gs_remain, tmp_data, tmp_sign2,
                                ind, with_weight, DEweights, tf_exprs, useCorSign)
        tmp_activity2 = tmp_rslt$activity
        tmp_sign2 = tmp_rslt$sign
        cor1 = apply(tmp_data, 1, function(x) {
            tmpc = cor(x, tmp_activity, method = "spearman")
            return(tmpc)
        })
        cor2 = apply(tmp_data, 1, function(x) {
            tmpc = cor(x, tmp_activity2, method = "spearman")
            return(tmpc)
        })
        if (if_module == T){
            all_list[[tf]] = tmp_sign
            all_activity = rbind(all_activity, tmp_activity)
        }else{
            mean1 = mean(abs(cor1))
            mean2 = mean(abs(cor2))
            if (mean1 < mean2) {
                tmp_sign = tmp_sign2
                tmp_activity = tmp_activity2
            }
            all_list[[tf]] = tmp_sign
            all_activity = rbind(all_activity, tmp_activity)
        }
    }
    dimnames(all_activity) = list(tfs, colnames(data))
    return(list(all_list = all_list, all_activities = all_activity[complete.cases(all_activity), ]))
}

#' The core function to compute the activity profile of an TF
#' @param gs_remain a vector of target genes after filtering
#' @param tmp_data gene expression of target genes
#' @param tmp_sign sign of target genes (+1 for one group, -1 for the other)
#' @param ind Hill coefficient parameter used in the weighting factors (default: 1/5, recommend to use 0 < ind < 1/4)
#' @param with_weight whether weighting factors (based on DEG p-values) are used to compute TF activity (default: TRUE)
#' @param DE_weights a vector of the input for computing DE weighting factors (typically, adjusted p-values from DEG analysis)
#' @param tf_exprs a vector of gene expression of the TF
#' @param useCorSign allow the option to use the TF gene expression correlation to flip signs (default: TRUE)
#' @return a list of results: 
#'         activity: matrix of TF activity.
#'         sign: grouping scheme of all TF gene sets.
cal_activity = function (gs_remain, tmp_data, tmp_sign, ind, with_weight, DE_weights, tf_exprs, useCorSign = useCorSign) {
    if(length(gs_remain) == 1) {
        tmp_activity = as.vector(scale(tmp_data[gs_remain,]))  # by default it's a matrix
        tmp_sign = tmp_sign[gs_remain, ,drop = FALSE ]
    } else{
        tmp_rem_data = row_norm(tmp_data[gs_remain, ])
        tmp_sign = tmp_sign[gs_remain, , drop = FALSE]
        if(with_weight) {
            tmp_weight = Hill(DE_weights[gs_remain, ], ind)
        } else {
            tmp_weight = rep(1.0, length(gs_remain))
        }
        tmp_activity = apply(tmp_rem_data, 2, function(x){
            rlst = x * tmp_sign * tmp_weight
            return(sum(rlst)/sum(tmp_weight))})
        tmp_activity = as.numeric(tmp_activity)
    }
  if(useCorSign){
    if (cor(tmp_activity, tf_exprs, method="spearman") < 0){
        tmp_activity = -tmp_activity
        tmp_sign = -tmp_sign
      # dimnames(tmp_sign) = gs_remain
    }
  }
    else{
        tmp_activity = tmp_activity
        tmp_sign = tmp_sign
      # dimnames(tmp_sign) = gs_remain
    }
    return(list(activity = tmp_activity, sign = tmp_sign))
}