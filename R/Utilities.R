####################################################
################ Utlities Functions ################
####################################################

#1. Activities Calculation
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

#2. remove Non-informative genes
rem_data = function(x){
    rem_ids = which(apply(x, 1, sd) == 0)
    if (length(rem_ids) > 0){
        return(x[-rem_ids,])
    }
    else{
        return(x)
    }
}

#3. Row normalization
row_norm = function(data){
    row_mean = apply(data, 1, mean)
    row_sd = apply(data, 1, sd)
    norm_data = apply(data, 2, function(x) (x - row_mean)/(row_sd))
    return(norm_data)
}

#4. Hill function for the gene weight
Hill = function(x, ind){
    return(1/ (1 + (x/0.05) ^(ind)))
}


#5. Downstream Grouping Results
Down_Streams = function(gs_str, gs_db, eset, cellwidth = 4, cellheight = 4, showname = T, use_module=T){
    require(pheatmap)
    require(Biobase)
    require(Hmisc)
    data = exprs(eset)
    g_inds = which(gs_db[[gs_str]] %in% rownames(data))
    tmp_data = data[gs_db[[gs_str]][g_inds], ]; rownames(tmp_data) = gs_db[[gs_str]][g_inds]
    cor_mat = rcorr(as.matrix(t(tmp_data)),type = "pearson")[['r']]; cor_mat = data.frame(cor_mat)
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

Down_Streams_nofiltering = function(gs_str, gs_db, eset, cellwidth = 4, cellheight = 4, showname = T){
    require(pheatmap)
    require(Biobase)
    require(Hmisc)
    if (is(eset, "ExpressionSet")){
        data = exprs(eset)
    }else{data = eset}
    g_inds = which(gs_db[[gs_str]] %in% rownames(data))
    tmp_data = data[gs_db[[gs_str]][g_inds], ]; rownames(tmp_data) = gs_db[[gs_str]][g_inds]
    cor_mat = rcorr(as.matrix(t(tmp_data)),type = "pearson")[['r']]; cor_mat = data.frame(cor_mat)

    pheatmap(cor_mat, show_colnames= showname ,show_rownames = showname,
             color = colorRampPalette(c("blue", "white", "red"))(20),
             main = gs_str, border_color = FALSE, cellwidth = cellwidth, cellheight = cellheight)
}

#7 Combine Activity and Expression Heatmap
Combine_heatmap = function(new_activity, eset){
    require(ComplexHeatmap)
    require(circlize)
    if (is(eset, "ExpressionSet")){
        data = exprs(eset)
    }else{data = eset}
    H1 = Heatmap(row_norm(new_activity), col = colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                 cluster_columns = F, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity")
    gs = rownames(new_activity)
    gc = colnames(new_activity)
    H2 = Heatmap(row_norm(data[gs, gc]), col = colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                 cluster_columns = F, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression")
    H1 + H2
}

# with Column clustering
Combine_heatmap2 = function(new_activity, eset){
    require(gplots)
    data = exprs(eset)
    dist.pear <- function(x) as.dist(1-cor(t(x)))
    hclust.ave <- function(x) hclust(x, method="average")

    my_palette <- colorRampPalette(c("green", "white", "red"))(n = 1000)
    colors = c(seq(-2, 2, length=1001))

    act = new_activity
    H1 = heatmap.2(row_norm(act), distfun =dist.pear, hclustfun =hclust.ave, col=my_palette, breaks = colors,
                   main = "Activity",  tracecol=NULL)
    gs = rownames(act)
    gc = colnames(act)

    exp = data[gs, gc]
    exp2 <- exp[rev(H1$rowInd),H1$colInd]
    H2 = heatmap.2(row_norm(as.matrix(exp2)), tracecol=NULL, col=my_palette, breaks = colors,
                   dendrogram = "none", Rowv = FALSE, Colv = FALSE, main = "Expression")
}

plot_network = function(tf_links){
    require(igraph)
    tf_nodes = unique(c(as.character(tf_links$from),
                        as.character(tf_links$to)))
    tf_graph = graph_from_data_frame(d=tf_links, vertices=tf_nodes, directed=T)
    cat("# of nodes: ", length(tf_nodes), " # of links: ", dim(tf_links)[1], "\n")
    #l = layout_with_fr(tf_graph)
    l = layout_with_lgl
    V(tf_graph)$size = 10
    types = floor((E(tf_graph)$relation+1)/2)   # 1 - activation; 2 - inhibition
    colors = E(tf_graph)$relation - types * 2 + 2   # 1 - transcription; 2 - degradation; 3 - signaling
    plot(tf_graph, edge.arrow.size= 0.2, edge.curved=0.2,
         edge.color= c("red", "blue")[colors],
         edge.lty = c("solid", "dashed","dotted")[types],
         vertex.color= "lightgreen",
         vertex.frame.color= NA,
         vertex.label=V(tf_graph)$from, vertex.label.color="black",
         vertex.label.cex=0.7,  vertex.label.font=1.0,
         layout = l)
}

## Convert to log10 (CPM) measurement in the RNA-Seq matrix
toCPM = function(ctMat){
    L = colSums(ctMat)/10e6
    cts = sweep(ctMat, 2, L, FUN="/")
    mat = log10(as.matrix(cts)+1)
    return(mat)
}

plot_network_v = function(tf_links = tf_links){
  require(visNetwork)
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
