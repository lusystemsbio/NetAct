########################################################
################ TF Activities Function ################
########################################################

TF_Activity = function (tfs, GSDB, eset, DErslt, with_weight = TRUE, if_module = FALSE,
                          ind = 1/5, useDatabaseSign = FALSE, useCorSign = TRUE, databaseFilename = "data/hDB_v2_sign_GTEx.RDS"){
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
    if(useDatabaseSign){
    signs <- readRDS(databaseFilename)
    print("Using Database sign")
    print(dim(signs))
    }
    for (tf in tfs) {
      # tf <- tfs[1]
       print(tf)
        comms = intersect(GSDB[[tf]], rownames(table))
        g_inds = which(rownames(data) %in% comms)
        tmp_data = data[g_inds, ]
        tmp_data = rem_data(tmp_data)
        # cor_mat = rcorr(as.matrix(t(tmp_data)), type = "spearman")[["r"]]
        cor_mat = cor(as.matrix(t(tmp_data)), method = "spearman")
        # tmp_dist <- cor_mat
        tmp_dist <- (cor_mat+1)/2
        diag(tmp_dist) <- 0 #(-(rowSums(mod_mat) -1))
        dimnames(cor_mat) = list(rownames(tmp_data), rownames(tmp_data))
        tmp_names = list(rownames(tmp_data), "sign")
        col_sums = apply(tmp_dist, 2, sum)
        mod_mat = tmp_dist - sapply(col_sums, function(x) x *
                                        col_sums/sum(tmp_dist))
        # mod_mat = cor_mat
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
        if(useDatabaseSign){
        #  print("Using Database sign")
          # print(head(signs))
          # print(tf)
          tfSign <- signs[which(signs$SrcGene == tf),]
          # print(tfSign)
          tfSign <- tfSign[tfSign$TgtGene %in% gs_remain,]
          negAvg <- mean(tfSign$Sign[which(tfSign$TgtGene %in% rownames(tmp_sign)[tmp_sign<0])])
          posAvg <- mean(tfSign$Sign[which(tfSign$TgtGene %in% rownames(tmp_sign)[tmp_sign>0])])
          if(!is.finite(posAvg)) posAvg <- 0
          if(!is.finite(negAvg)) negAvg <- 0
          # print(posAvg - negAvg)
          # print(posAvg)
          if(is.finite(negAvg) & is.finite(posAvg) & ((posAvg - negAvg)<0) ){ tmp_sign = -tmp_sign
          print(tf)
          }
        }
        
        tmp_rslt = cal_activity(gs_remain, tmp_data, tmp_sign,
                                ind, with_weight, DEweights, tf_exprs, useCorSign)
        #tmp_activity[which(kdGene=="EP300")]
        
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
        if(useDatabaseSign){
          tfSign <- signs[signs$SrcGene == tf,]
          tfSign <- tfSign[tfSign$TgtGene %in% gs_remain,]
          negAvg <- mean(tfSign$Sign[which(tfSign$TgtGene %in% rownames(tmp_sign)[tmp_sign<0])])
          posAvg <- mean(tfSign$Sign[which(tfSign$TgtGene %in% rownames(tmp_sign)[tmp_sign>0])])
          if(is.finite(negAvg) & is.finite(posAvg) & ((posAvg - negAvg)<0) ) tmp_sign = -tmp_sign
        }
        
       # write(gs_remain,paste0(tf,".txt"))
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
    return(list(all_list = all_list, all_activities = all_activity[complete.cases(all_activity),
                                                                   ]))
}
