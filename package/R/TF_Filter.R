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

## Filter the links in the topology; 1 means activation 2 means inhibition
TF_Filter = function(actMat, GSDB, miTh = 1.4, nbins = 8, method = "spearman"){
    require(infotheo)
    NetMat = allNet(GSDB)
    # calculate the MI
    miMat = discretize(t(actMat), disc="equalfreq", nbins = 8)
    miMat = mutinformation(miMat, method = "shrink")
    diag(miMat) = 0
    actLinks = melt(miMat)
    corMat = cor(t(actMat), method = method)
    corLinks = melt(corMat)
    
    # construct the links
    tf_source = actLinks[which(actLinks$value > miTh),1]
    tf_target = actLinks[which(actLinks$value > miTh),2]
    tf_links = data.frame(as.matrix(NetMat[NetMat$from %in% tf_source & NetMat$to %in% tf_target, ]), 
                          stringsAsFactors = F)
    tf_links$relation = NULL
    
    # determine the directions
    for (i in 1:nrow(tf_links)){
        source_name = tf_links$from[i]
        target_name = tf_links$to[i]
        tf_links$relation[i] <- ifelse(corMat[source_name, target_name] > 0, 1, 2)
    }
    rownames(tf_links) = NULL
    return(tf_links)
}

