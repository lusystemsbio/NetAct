####################################################
################ Utlities Functions ################
####################################################

#' Remove Non-informative genes
#' @param x: gene expression matrix
#' @return x: gene expression matrix without containing non-informative genes
rem_data = function(x){
    rem_ids = which(apply(x, 1, sd) == 0)
    if (length(rem_ids) > 0){
        return(x[-rem_ids,])
    }
    else{
        return(x)
    }
}

#' Row normalization (standardization)
#' @param data: gene expression matrix
#' @return norm_data: standardized gene expression matrix
row_norm = function(data){
    row_mean = apply(data, 1, mean)
    row_sd = apply(data, 1, sd)
    norm_data = apply(data, 2, function(x) (x - row_mean)/(row_sd))
    return(norm_data)
}

#' Hill function for the gene weight
#' @param x: value (adj p-value)
#' @param ind: Hill coefficient
#' @return Hill function of x
Hill = function(x, ind){
    return(1/ (1 + (x/0.05) ^(ind)))
}

#' Downstream Grouping Results
#' @param gs_str: Vector. A list of gene set names (TFs)
#' @param gs_db: List of list. Gene set data base
#' @param eset: ExpressionSet of gene expression data
#' @param cellwidth: plot width (default 4)
#' @param cellheight:  plot height (default 4)
#' @param showname: show colnames and rownames (default TRUE)
#' @param use_module: use modularity of gene grouping (default TRUE)
#' @return List: gene grouping scheme
#' @import pheatmap
#' @import Biobase
#' @export
Down_Streams = function(gs_str, gs_db, eset, cellwidth = 4, cellheight = 4, showname = T, use_module=T){
#    require(pheatmap)
#    require(Biobase)
#    require(Hmisc)
    data = exprs(eset)
    g_inds = which(gs_db[[gs_str]] %in% rownames(data))
    tmp_data = data[gs_db[[gs_str]][g_inds], ]; rownames(tmp_data) = gs_db[[gs_str]][g_inds]
#    cor_mat = rcorr(as.matrix(t(tmp_data)),type = "pearson")[['r']]; cor_mat = data.frame(cor_mat)
    cor_mat = cor(as.matrix(t(tmp_data)), method = "spearman")
    cor_mat = data.drame(cor_mat)
    dimnames(cor_mat) = list(rownames(tmp_data), rownames(tmp_data))
    tmp_names = list(rownames(tmp_data), "sign")

    # Use the correlation for the pairwise relation
    tmp_dist = cor_mat
    col_sums = apply(tmp_dist, 1, sum)
    mod_mat = tmp_dist - sapply(col_sums, function(x) x*col_sums/ (sum(tmp_dist)))
    ## Find the largest eigen value and corresponding vector
    ev = eigen(mod_mat)
    tmp_sign = data.frame(sign(ev$vectors[,1])); dimnames(tmp_sign) = tmp_names
    ng_mat = cor_mat[tmp_sign$sign == -1, tmp_sign$sign == -1]
    ps_mat = cor_mat[tmp_sign$sign == 1, tmp_sign$sign == 1]

    if (use_module ==F){
        gs_list = list(); gs_remain = c()
        tf_exprs = data[tf, ]
        cors = apply(tmp_data, 1, function(x) cor(x, tf_exprs,
                                                  method = "spearman"))
        pos_cors = cors[cors > 0]
        neg_cors = cors[cors < 0]
        tmp_sign2 = data.frame(sign(cors))
        m_pos = mean(pos_cors)
        m_neg = mean(neg_cors)
        gs_list[[1]] = names(pos_cors[pos_cors >= m_pos])
        gs_list[[2]] = names(neg_cors[neg_cors <= m_neg])
        gs_remain = unlist(gs_list)
        tmp_rslt = cal_activity(gs_remain, tmp_data, tmp_sign2,
                                ind, with_weight, DE_weights, tf_exprs)
        tmp_activity2 = tmp_rslt$activity

    }else{
        # Use Kmeans to filter outlies within one module
        test_mats = list(ng_mat = ng_mat, ps_mat = ps_mat)
        gs_remain = c(); gs_list = list()
        for (i in 1:2){
            test_mat = test_mats[[i]]
            if (nrow(test_mat) > 1){
                tmp_kmean = kmeans(test_mat, 1)
                centers = tmp_kmean$centers
                # use correlation distance
                distances = apply(test_mat, 1, function(x) cor(as.numeric(x), as.numeric(centers),method = "pearson"))
                m = mean(distances[distances > 0])
                gs_list[[i]] = names(distances)[distances > m]
                gs_remain = c(gs_remain, names(distances)[distances > m])
            }
        }
        if(length(gs_list) == 1){
            pheatmap(cor_mat[gs_remain, gs_remain], show_colnames= F ,show_rownames = F,
                     color = colorRampPalette(c("blue", "white", "red"))(20),
                     main = gs_str, border_color = FALSE, cellwidth = cellwidth, cellheight = cellheight)
        }
        else{
            annotation_rc = data.frame(genes = c(rep("Type1", length(gs_list[[1]])),
                                                 rep("Type2", length(gs_list[[2]]))))
            rownames(annotation_rc) = gs_remain
            genes =  c("purple3", "orange3"); names(genes) <- c("Type1", "Type2")
            anno_colors = list(genes = genes)
            pheatmap(cor_mat[gs_remain, gs_remain], show_colnames= showname ,show_rownames = showname,
                     annotation_row = annotation_rc, annotation_col = annotation_rc,
                     annotation_colors = anno_colors,
                     color = colorRampPalette(c("blue", "white", "red"))(20),
                     main = gs_str, border_color = FALSE, cellwidth = cellwidth, cellheight = cellheight)
        }
    }
    return(gs_list)
}

#' Downstream Grouping Results, without doing gene filtering
#' @param gs_str: Vector. A list of gene set names (TFs)
#' @param gs_db: List of list. Gene set data base
#' @param eset: ExpressionSet of gene expression data
#' @param cellwidth: plot width (default 4)
#' @param cellheight:  plot height (default 4)
#' @param showname: show colnames and rownames (default TRUE)
#' @return List: gene grouping scheme
#' @import pheatmap
#' @import Biobase
#' @export
Down_Streams_nofiltering = function(gs_str, gs_db, eset, cellwidth = 4, cellheight = 4, showname = T){
#    require(pheatmap)
#    require(Biobase)
#    require(Hmisc)
    if (is(eset, "ExpressionSet")){
        data = exprs(eset)
    }else{data = eset}
    g_inds = which(gs_db[[gs_str]] %in% rownames(data))
    tmp_data = data[gs_db[[gs_str]][g_inds], ]; rownames(tmp_data) = gs_db[[gs_str]][g_inds]
#    cor_mat = rcorr(as.matrix(t(tmp_data)),type = "pearson")[['r']]; cor_mat = data.frame(cor_mat)
    cor_mat = cor(as.matrix(t(tmp_data)), method = "spearman")
    cor_mat = data.drame(cor_mat)

    pheatmap(cor_mat, show_colnames= showname ,show_rownames = showname,
             color = colorRampPalette(c("blue", "white", "red"))(20),
             main = gs_str, border_color = FALSE, cellwidth = cellwidth, cellheight = cellheight)
}

#' Plotting TF gene expresion & activity heatmap
#' @param new_activity: Matrix. TF activity matrix
#' @param eset: ExpressionSet of gene expression data
#' @return Heatmap plotting object
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @export
Activity_heatmap = function(new_activity, eset){
#    require(ComplexHeatmap)
#    require(circlize)
    if (is(eset, "ExpressionSet")){
        data = exprs(eset)
    }else{data = eset}
    H1 = Heatmap(row_norm(new_activity), col = colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                 cluster_columns = T, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity")
    gs = rownames(new_activity)
    gc = colnames(new_activity)
    H2 = Heatmap(row_norm(data[gs, gc]), col = colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                 cluster_columns = T, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression")
    H1 + H2
}

#' convert to log10 (CPM) measurement in the RNA-Seq matrix
#' @param ctMat: Matrix of gene expression counts
#' @return mat: Matrix of CPM gene expression
#' @export
toCPM = function(ctMat){
    L = colSums(ctMat)/10e6
    cts = sweep(ctMat, 2, L, FUN="/")
    mat = log10(as.matrix(cts)+1)
    return(mat)
}

#' Plotting gene network
#' @param tf_links: a data frame of networ interactions
#' @return visNetwork object
#' @import visNetwork
#' @export
plot_network = function(tf_links = tf_links){
#  require(visNetwork)
  topology=data.frame(as.matrix(tf_links), stringsAsFactors = F)

  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list, label = node_list, font.size =30, value=c(rep(1,length(node_list))))

  #nodes <- data.frame(id = node_list, label = node_list, font.size =30,shape='circle',value=c(rep(1,length(node_list))))
  edge_col <- data.frame(c(1,2),c("blue","darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1,2),c("arrow","circle"))
  colnames(arrow_type) <- c("type", "color")
  edges <- data.frame(from =c(topology[,1]), to = c(topology[,2])
                      , arrows.to.type	=arrow_type$color[c(as.numeric(topology[,3]))]
                      , width = 3
                      , color = edge_col$color[c(as.numeric(topology[,3]))]
  )
  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
#  file  <- paste("network_",file,".html",sep="")
#  visSave(network, file = file, selfcontained = F)
}

#' Filtered gene set database based on minimum sizes
#' @param GSDB: list of list. gene set database
#' @param geneList: a vector of available genes
#' @param minSize: minimum number of genes of a gene set (default: 5)
#' @return DB: list of list. filtered gene set database
#' @export
filterDB <- function(GSDB,geneList,minSize = 5)
{
  nameTF <- names(GSDB)
  nameTF <- nameTF[nameTF %in% geneList]
  nTF <- length(nameTF)
  nameTF <- names(DB)
  j=1
  for(i in 1:nTF){
    if(nameTF[i] %in% geneList){
      DB[[j]] <- DB[[j]][DB[[j]] %in% geneList]
      if(length(DB[[j]]) > (minSize-1))
        {
      j=j+1
      }
      else{(DB[j] <- NULL)}
    }
    else{
      DB[j] <- NULL
      
    }
  }
  return(DB)
}
