###############################################################################################
#################################### Filter links Function #################################### 
###############################################################################################

## Compile a GSDB to a matrix with 2 columns
allNet = function(GSDB){
    allNet = matrix(unlist(GSDB), ncol = 1)
    allNet = cbind(rep(names(GSDB), as.vector(sapply(GSDB, length))), allNet)
    allNet = data.frame(allNet); colnames(allNet) = c("from", "to")
    return(allNet)
}

#' @export
#' @title Calculate mutual information
#' @description Mutual information between all pairs based on entropy package.
#' @param actMat numeric matrix. 
#' @param nbins integer (optional) Number of bins Default 16
#' @return numeric matrix (0-1) Matrix containing mutual information values
#' 
calculateMI <- function(actMat = actMat, nbins=16){
  nGenes <- dim(actMat)[1]
  miMat <- matrix(0,nrow = nGenes,ncol = nGenes)
  geneNames <- rownames(actMat)
  rownames(miMat) <- geneNames
  colnames(miMat) <- geneNames
  
  #i=1
  for(i in 1:(nGenes-1))
  {
    for(j in (i+1):(nGenes))
    {
      temp <- entropy::discretize2d(actMat[i,],actMat[j,],nbins,nbins)
      miMat[i, j] <- entropy::mi.shrink(temp,verbose = F)
      miMat[j, i] <- miMat[i, j]
    }
  }
  return(miMat)
}


## Filter the links in the topology; 1 means activation 2 means inhibition
#' @export
#' @title Generate network
#' @description Network calculated using activity and interaction database. Uses
#' mutual information to find possible interactions and keeps the interactions 
#' if they are available in the database. Sign of interaction is assigned based
#' on the correlation between the activities.
#' @param actMat numeric matrix. Matrix containing the activities
#' @param GSDB List of list Gene set database of interactions
#' @param maxTf integer (optional) Default 75. Maximum number of transcription 
#' factors in the network. If \code{removeSignalling} is \code{TRUE} 
#' the actual number will be less. 
#' @param maxInteractions integer (optional) Default 300. Maximum number of
#' interactions in the network.
#' @param miTh numeric. Mutual information threhold 
#' @param nbins integer (optional) Number of bins Default 16.
#' @param corMethod character (optional) Method to compute correlaiton. See 
#' \link{\code[stats]{cor}}.
#' @param useCor Logical (optional) Whether to use correlation instead of
#' mutual information to find possible interactions. Default \code{FALSE}.
#' @param removeSignalling logical (optional). Whether to remove the Tfs which
#' are not the target of any other Tfs. Default TRUE. It is not recursive and 
#' the generated network might still contain some signallings tfs.
#' @param DPI logical (optional) Default TRUE. 
#' Whether to apply the data processing
#' inequality to remove weak edges from triangles. 
#' @return data.frame Contains the interactions in a dataframe listing 
#' source tf, target tf and interaction type (1-activation, 2-inhibition).  
#' 
#' 
TF_Filter = function(actMat, GSDB, miTh = 0.4, maxTf = 75, 
                     maxInteractions = 300,  
                     nbins = 16, corMethod = "spearman", useCor = FALSE, 
                     removeSignalling = FALSE, DPI = FALSE, ...){

  corMat = cor(t(actMat), method = corMethod)
  if(useCor){
    miMat <- abs(corMat)
    diag(miMat) <- 0
  }
  else{
    miMat <- calculateMI(actMat, nbins)
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
    
    return(tfLinks)

}

#' @export
#' @title Apply data processing inequality
#' @description Remove the interactions from a triangle which have lowest 
#' interaction score.
#' @param tfLinks Dataframe containing the interactions as source (character),
#'  target (character), type (integer).
#' @param miMat numeric matrix Interaction scores based on mutual information or 
#' correlation. 
#' @param miDiff numeric (0-1) Default 0.0 (optional) Minimum difference 
#' between mutual informations of a traingle for the edge to be removed.
#' @param minMiTh numeric (0-1) Default 0.5. Minimum value of MI for an interaction 
#' which will not be removed. 
#' @return data.frame containing the filtered interactions.
#'   
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