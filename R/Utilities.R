####################################################
################ Utlities Functions ################
####################################################

#' Remove Non-informative genes
#' @param x gene expression matrix
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
#' @param data gene expression matrix
#' @return norm_data: standardized gene expression matrix
row_norm = function(data){
    row_mean = apply(data, 1, mean)
    row_sd = apply(data, 1, sd)
    norm_data = apply(data, 2, function(x) (x - row_mean)/(row_sd))
    return(norm_data)
}

#' Hill function for the gene weight
#' @param x value (adj p-value)
#' @param ind Hill coefficient
#' @return Hill function of x
Hill = function(x, ind){
    return(1/ (1 + (x/0.05) ^(ind)))
}

#' Plotting TF gene expresion & activity heatmap
#' @param new_activity Matrix. TF activity matrix
#' @param eset ExpressionSet of gene expression data
#' @return Heatmap plotting object
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap 
#' @export
Activity_heatmap = function(new_activity, eset){
#    require(circlize)
    if (is(eset, "ExpressionSet")){
        data = exprs(eset)
    }else{data = eset}
    H1 = Heatmap(row_norm(new_activity), col = colorRamp2(c(-2, 0, 2), c("blue3", "white", "red")),
                 cluster_columns = T, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity")
    gs = rownames(new_activity)
    gc = colnames(new_activity)
    H2 = Heatmap(row_norm(data[gs, gc]), col = colorRamp2(c(-2, 0, 2), c("blue3", "white", "red")),
                 cluster_columns = T, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression")
    H1 + H2
}

#' convert to log10 (CPM) measurement in the RNA-Seq matrix
#' @param ctMat Matrix of gene expression counts
#' @return mat: Matrix of CPM gene expression
#' @export
toCPM = function(ctMat){
    L = colSums(ctMat)/10e6
    cts = sweep(ctMat, 2, L, FUN="/")
    mat = log10(as.matrix(cts)+1)
    return(mat)
}

#' Plotting gene network
#' @param tf_links a data frame of networ interactions
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
#' @param GSDB list of list. gene set database
#' @param geneList a vector of available genes
#' @param minSize minimum number of genes of a gene set (default: 5)
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
