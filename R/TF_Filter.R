###############################################################################################
#################################### Filter links Function #################################### 
###############################################################################################

#' Compile a GSDB to a matrix with 2 columns
#' @param GSDB List of list. Gene set database of interactions
#' @return matrix. Matrix, each row containing regulators ("from"), targets ("to")
#' @export
allNet = function(GSDB){
    allNet = matrix(unlist(GSDB), ncol = 1)
    allNet = cbind(rep(names(GSDB), as.vector(sapply(GSDB, length))), allNet)
    allNet = data.frame(allNet); colnames(allNet) = c("from", "to")
    return(allNet)
}

#' @title Calculate mutual information
#' @description Mutual information between all pairs based on entropy package.
#' @param actMat numeric matrix. 
#' @param nbins integer (optional). Number of bins Default 16
#' @param method MI calculation method: 
#'                "2d": 2D discretization with entropy (default) 
#'                "1d": 1D discretization with infotheo
#' @import entropy
#' @import infotheo
#' @return numeric matrix (0-1). Matrix containing mutual information values
calculateMI <- function(actMat = actMat, nbins=16, method = "2d"){
  nGenes <- dim(actMat)[1]
  geneNames <- rownames(actMat)
  
#  mat_cor = cor(t(actMat))
  
  if(method == "1d"){
    temp = infotheo::discretize(t(actMat), disc="equalfreq", nbins = nbins)
    miMat = infotheo::mutinformation(temp, method = "shrink") 
    rownames(miMat) <- geneNames
    colnames(miMat) <- geneNames
  }else{   # "2d", default choice
    miMat <- matrix(0,nrow = nGenes,ncol = nGenes)
    rownames(miMat) <- geneNames
    colnames(miMat) <- geneNames
    for(i in 1:(nGenes-1))
    {
      for(j in (i+1):(nGenes))
      {
        temp <- entropy::discretize2d(actMat[i,],actMat[j,],nbins,nbins)
        miMat[i, j] <- entropy::mi.shrink(temp,verbose = F)
        miMat[j, i] <- miMat[i, j]
      }
    }
  }
  
  #i=1
  return(miMat)
}

## Filter the links in the topology; 1 means activation 2 means inhibition
#' @title Generate network
#' @description Network calculated using activity and interaction database. Uses
#' mutual information to find possible interactions and keeps the interactions 
#' if they are available in the database. Sign of interaction is assigned based
#' on the correlation between the activities.
#' @param actMat numeric matrix. Matrix containing the activities
#' @param GSDB List of list. Gene set database of interactions
#' @param miTh numeric. Mutual information threshold 
#' @param maxTf integer (optional). Default 75. Maximum number of transcription 
#' factors in the network. If \code{removeSignalling} is \code{TRUE} 
#' the actual number will be less. 
#' @param maxInteractions integer (optional). Default 300. Maximum number of
#' interactions in the network.
#' @param nbins integer (optional). Number of bins Default 16.
#' @param miMethod MI calculation method: 
#'                "2d": 2D discretization with entropy (default) 
#'                "1d": 1D discretization with infotheo
#' @param corMethod character (optional). Method to compute correlation.
#' @param useCor Logical (optional). Whether to use correlation instead of
#' mutual information to find possible interactions. 
#' @param removeSignalling logical (optional). Whether to remove the Tfs which
#' are not the target of any other Tfs. Default TRUE. It is not recursive and 
#' the generated network might still contain some signalling tfs.
#' @param DPI logical (optional). Default FALSE. 
#' Whether to apply the data processing
#' inequality to remove weak edges from triangles. 
#' @param nameFile character (optional). Ouput file name. Default NULL (no file output). 
#' @param ... two additional parameters passed from applyDPI (default: miDiff = 0, minMiTh = 0.5)
#' @return data.frame. Contains the interactions in a dataframe listing.
#' source tf, target tf and interaction type (1-activation, 2-inhibition).  
#' @import reshape2
#' @export
TF_Filter = function(actMat, GSDB, miTh = 0.4, maxTf = 75, 
                     maxInteractions = 300, nbins = 16, miMethod = "2d",
                     corMethod = "spearman", useCor = FALSE, 
                     removeSignalling = FALSE, DPI = FALSE, nameFile = NULL, ...){
  
  corMat = cor(t(actMat), method = corMethod)
  if(useCor){
    miMat <- abs(corMat)
    diag(miMat) <- 0
  }
  else{
    miMat <- calculateMI(actMat, nbins, method = miMethod)
  }
  
  actLinks = reshape2::melt(miMat)
  actLinks <- actLinks[order(actLinks$value,decreasing = TRUE),]
  actLinks <- actLinks[actLinks$value > miTh,]
  actLinks <- data.frame(actLinks, stringsAsFactors = FALSE)
  tf <- 0
  link <- 0
  
  tfLinks <- as.data.frame(matrix(NA, ncol = 2,nrow = maxInteractions))
  names(tfLinks) <- c("Source", "Target")
  DBtf <- names(GSDB)
  counter <- 1
  while((link <= maxInteractions) && (tf <= maxTf) && 
        counter <= dim(actLinks)[1]){
    if(actLinks[counter,2] %in% GSDB[[which(DBtf == actLinks[counter,1])]]){
      link=link+1
      tfLinks[link, 1] <- as.character(actLinks[counter,1])
      tfLinks[link, 2] <- as.character(actLinks[counter,2])
      tf = length(union(tfLinks$Source,tfLinks$Target))
      
    }
    counter=counter+1
    
  }
  if(removeSignalling) {
    tfLinks <- tfLinks[tfLinks[,1] %in% tfLinks[,2],]
  }
  tfLinks <- tfLinks[complete.cases(tfLinks),]
  corLinks = reshape2::melt(corMat)
  
  tfLinks$Interaction = NULL
  
  # determine the directions
  for (i in 1:nrow(tfLinks)){
    tfLinks$Interaction[i] <- ifelse(corMat[tfLinks[i,1], 
                                            tfLinks[i,2]] > 0, 1, 2)
  }
  rownames(tfLinks) = NULL
  
  if(DPI){
    tfLinks <- applyDPI(tfLinks, miMat,...)
  }
  
  if (!is.null(nameFile)) {
    write.csv(tfLinks, file = paste0(nameFile, ".csv" ))
    saveRDS(tfLinks, file = paste0(nameFile, ".tpo" ))
  }
   
  if(nrow(tfLinks) == maxInteractions){
    print("Maximum number of interactions is reached!")
    
  }
  return(tfLinks)
}

#' @title Apply data processing inequality
#' @description Remove the interactions from a triangle which have lowest 
#' interaction score.
#' @param tfLinks Data.frame. containing the interactions as source (character),
#'  target (character), type (integer).
#' @param miMat numeric matrix. Interaction scores based on mutual information or 
#' correlation. 
#' @param miDiff numeric (0-1). Default 0.0 (optional) Minimum difference 
#' between mutual information of a traingle for the edge to be removed.
#' @param minMiTh numeric (0-1). Default 0.5. Minimum value of MI for an interaction 
#' which will not be removed. 
#' @return data.frame. containing the filtered interactions.
applyDPI <- function(tfLinks = tfLinks, miMat = miMat, miDiff = 0, minMiTh = 0.5){
  print(miDiff)
  print(minMiTh)
  edgeKeep <- vector(mode = "integer", length = length(tfLinks[,1]))
  edgeRemove <- vector(mode = "integer", length = length(tfLinks[,1]))
  for(i in seq_along(tfLinks[,1])){
    srcGene <- tfLinks[i,1]
    tgtGene <- tfLinks[i,2]
    allTgts <- tfLinks[which(tfLinks[,1] == srcGene),2]
    for(j in seq_along(allTgts)){
      threeGenes <- c(srcGene, tgtGene, allTgts[j])
      if(length(unique(threeGenes)) == 3){
        thirdGenes <- union(tfLinks[which(tfLinks[,1] == allTgts[j]),2],
                            tfLinks[which(tfLinks[,2] == allTgts[j]),1])
        if((srcGene %in% thirdGenes) & (tgtGene %in% thirdGenes)){
          tmp <- miMat[c(srcGene,tgtGene,allTgts[j]),c(srcGene, 
                                                       tgtGene,allTgts[j])]
          minValue <- which.min(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          minMi <- min(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          maxMi <- max(c(tmp[2,1],tmp[3,1],tmp[3,2]))
          edge <- list()
          edge[[1]] <- which((tfLinks[,1] == srcGene) & (tfLinks[,2] == tgtGene) |
                             (tfLinks[,2] == srcGene) & (tfLinks[,1] == tgtGene)
                           )
          edge[[2]] <- which((tfLinks[,1] == srcGene) & 
                             (tfLinks[,2] == allTgts[j]) | 
                             (tfLinks[,2] == srcGene) & 
                             (tfLinks[,1] == allTgts[j]))
          edge[[3]] <- which((tfLinks[,1] == tgtGene) & 
                             (tfLinks[,2] == allTgts[j]) | 
                             (tfLinks[,2] == tgtGene) & 
                             (tfLinks[,1] == allTgts[j]))

       #   edgeKeep[edge[1]] <- edgeKeep[edge[1]] + 1
      #    edgeKeep[edge[2]] <- edgeKeep[edge[2]] + 1
      #    edgeKeep[edge[3]] <- edgeKeep[edge[3]] + 1
          
       #   edgeKeep[edge[minValue]] <- edgeKeep[edge[minValue]] - 1
       #   edgeRemove[edge[minValue]] <- edgeRemove[edge[minValue]] - 1
          if((maxMi - minMi) > miDiff){
            if(minMi < minMiTh){
          edgeRemove[edge[[minValue]]] <- -1
            }
          }
        }
      }
    }
  }

  tfLinks <- tfLinks[-(which(edgeRemove < 0)),]
  return(tfLinks)
}

#' Obtain the adjacency matrix from a matrix of tf-target relationships
#' @param tfLinks matrix. Matrix of tf-target relationships
#' @return adjMat, matrix. adjacency matrix
getAdjacencyMat <- function(tfLinks = tfLinks){
  networkGenes <-  unique(c(tfLinks[,1],tfLinks[,2]))
  nGenes <- length(networkGenes)
  adjMat <- data.frame(matrix(data = 0,nrow = nGenes, ncol = nGenes))
  rownames(adjMat) <- networkGenes
  colnames(adjMat) <- as.character(networkGenes)
  for(i in 1:dim(tfLinks)[1]){
    adjMat[tfLinks[i,2], tfLinks[i,1]] <- tfLinks[i,3]
  }
  return(adjMat)
}

#' Generate network (an extension of TF_Filter)
#' @description Network calculated using activity and interaction database. Uses
#' mutual information to find possible interactions and keeps the interactions 
#' if they are available in the database. Sign of interaction is assigned based
#' on the correlation between the activities. An extension of TF_Filter. 
#' Add a list of genes of interest.
#' @param actMat numeric matrix. Matrix containing the activities
#' @param GSDB List of list. Gene set database of interactions
#' @param genes vector. a vector of gene symbols of genes of interest
#' @param DEgenes vector. a vector of gene symbols of DE genes 
#' @param eset expression set of gene expression data or gene expression matrix
#' @param miTh numeric. Mutual information threshold 
#' @param maxTf integer (optional). Default 75. Maximum number of transcription 
#' factors in the network. If \code{removeSignalling} is \code{TRUE} 
#' the actual number will be less. 
#' @param maxInteractions integer (optional). Default 300. Maximum number of
#' interactions in the network.
#' @param nbins integer (optional). Number of bins Default 16.
#' @param miMethod MI calculation method: e: entropy (default) or i: infotheo
#' @param corMethod character (optional). Method to compute correlation.
#' @param useCor Logical (optional). Whether to use correlation instead of
#' mutual information to find possible interactions. Default FALSE 
#' @param removeSignalling logical (optional). Whether to remove the Tfs which
#' are not the target of any other Tfs. Default FALSE. It is not recursive and 
#' the generated network might still contain some signalling tfs.
#' @param DPI logical (optional). Default FALSE. 
#' Whether to apply the data processing
#' inequality to remove weak edges from triangles. 
#' @param ... two additional parameters passed from applyDPI (default: miDiff = 0, minMiTh = 0.5)
#' @return List of data.frame. Contains the interactions in a data frame listing.
#' source tf, target tf and interaction type (1-activation, 2-inhibition).  
#'        tf_links: network interactions.
#'        new_links: new interactions associated with the genes of interest.
#' @export
TF_Filter_addgene <- function(actMat, GSDB, genes, DEgenes, eset, miTh = 0.4, maxTf = 75, 
                       maxInteractions = 300, nbins = 16, miMethod = "e",
                       corMethod = "spearman", useCor = FALSE, 
                       removeSignalling = FALSE, DPI = FALSE, ...){
    
# genes : the genes we want check 
#  DE genes: diffentially expressed gene list(get from the function "MicroDegs")
#  eset: the whole data set with all genes 
  if (is(eset, "ExpressionSet")){
    data = exprs(eset)
  }else{
    data = eset
  }
  genelist=intersect(genes,DEgenes)
  genelist=setdiff(genelist,names(GSDB))
  
  if(length(genelist) > 0){
    GSDB2=list()
    for (i in 1: length(genelist)){
      GSDB2[[i]]=c("NULL")
    }
    names(GSDB2)=genelist
    GSDB.n=append(GSDB,GSDB2)
  }else{
    print("No new genes added!")
    return(list(tf_links=tf_links2, new_links = NULL))
  }
  
  k1=as.matrix(rbind(actMat,data[genelist,]))
  tf_links = TF_Filter(actMat, GSDB, miTh =miTh, maxTf =maxTf, maxInteractions=maxInteractions, nbins =nbins, miMethod = miMethod, 
                       corMethod =corMethod, useCor=useCor,removeSignalling=removeSignalling,  DPI = DPI, nameFile = NULL,...)
  tf_links2 = TF_Filter(k1, GSDB.n, miTh =miTh, maxTf =maxTf, maxInteractions=maxInteractions, nbins =nbins, miMethod = miMethod,
                       corMethod =corMethod, useCor=useCor,removeSignalling=removeSignalling,  DPI = DPI, nameFile = NULL,...)
  new=setdiff(tf_links2[,2],tf_links[,2])
  tf_link_new=tf_links2[tf_links2[,2]%in%new,]
  return(list(tf_links=tf_links2,new_links=tf_link_new))
}