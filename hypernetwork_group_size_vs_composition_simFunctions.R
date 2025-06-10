##Generate individuals in family/kinship groups and spatially embed them in x-y grid
generate_population <- function(popSize, nInitInformed) {
  #Assign IDs
  id <- seq(1, popSize, 1)
  
  #Assign individual coordinates
  ind_x <- ind_y <- rep(NA, popSize)
  for(i in 1:popSize) {
    ind_x[i] <- runif(1, 0, 1)
    ind_y[i] <- runif(1, 0, 1)
  }
  
  ind_data <- data.frame(id,ind_x,ind_y)
  
  #Create an initial informed individual
  ind_data$informed <- 0
  ind_data$informed[sample(1:nrow(ind_data), nInitInformed, replace = FALSE)] <- 1
  
  #Identify initial informed nodes so that initial seed can be denoted
  informedNodes <- ind_data[which(ind_data$informed == 1),]$id
  ind_data$initInformed <- 0
  ind_data$initInformed[informedNodes] <- 1
  ind_data$acqTime <- 0
  ind_data$firstProd <- 0
  
  return(ind_data)
}

##Create spatial proximity network and hypergraphs
generate_latent_space_multilayer_hypergraph <- function(ind_data, r) {
  #Calculate distance matrix
  dist_mat <- as.matrix(dist(cbind(ind_data$ind_x, ind_data$ind_y)))
  diag(dist_mat) <- 0
  
  #Set initial threshold to use for calculating dyadic layer and intermediate layer (i.e., hypernetwork 1)
  thresh <- r[1]
  
  #Create matrix of dyadic associations
  sp_mat <- dist_mat < thresh
  sp_net <- apply(sp_mat, 2, as.numeric)
  diag(sp_net) <- 0
  sp_graph <- graph_from_adjacency_matrix(sp_net, mode="max")
  
  #The following procedure ensures that a single connected component is formed
  C1 <- components(sp_graph)
  numComps <- max(C1$membership)
  
  if(numComps > 1) {
    
    #Create data.table to collect nearest distances among network components
    newEdgeList <- data.table("Comp1" = unlist(sapply(1:(numComps-1), function(x) rep(x,(numComps-x)))),
                              "Comp2" = unlist(sapply(1:(numComps-1), function(x) (x+1):numComps)),
                              "ID1" = 0,
                              "ID2" = 0,
                              "Dist" = 0)
    
    #For each pair of network components, identify the pair of nodes that are physically closest
    for(i in 1:(numComps-1)) {
      for(j in 2:numComps) {
        memb1 <- c(which(C1$membership == i))
        memb2 <- c(which(C1$membership == j))
        tempMat <- sapply(memb1, function(x) dist_mat[memb2,x])
        newEdgeList[Comp1 == i & Comp2 == j]$Dist <- min(tempMat)
        targetEdge <- which(tempMat == min(tempMat))
        if(C1$csize[i] > 1 & C1$csize[j] > 1){
          newEdgeList[Comp1 == i & Comp2 == j]$ID1 <- as.vector(memb1[ceiling(targetEdge/dim(tempMat)[1])])
          newEdgeList[Comp1 == i & Comp2 == j]$ID2 <- as.vector(memb2[targetEdge-(dim(tempMat)[1] * (ceiling(targetEdge/dim(tempMat)[1])-1))])
        } else{
          if(C1$csize[i] > 1 & C1$csize[j] == 1) {
            newEdgeList[Comp1 == i & Comp2 == j]$ID2 <- memb2
            newEdgeList[Comp1 == i & Comp2 == j]$ID1 <- memb1[targetEdge]
          } else {
            newEdgeList[Comp1 == i & Comp2 == j]$ID1 <- memb1
            newEdgeList[Comp1 == i & Comp2 == j]$ID2 <- memb2[targetEdge]
          }
        }
      }
    }
    newEdgeList <- newEdgeList[order(newEdgeList$Dist)]
    
    #Iterate through the new potential edges in order of shortest distance
    #A new edge is added one at a time and the number of components in the resulting network assessed
    #Stop the process once a single connected component has been formed
    i = 1
    repeat{
      sp_net[newEdgeList$ID1[i], newEdgeList$ID2[i]] <- sp_net[newEdgeList$ID2[i], newEdgeList$ID1[i]] <- 1
      if(components(graph_from_adjacency_matrix(sp_net, mode ="max"))$no == 1) {
        break}
      i <- i + 1
    }
    
    #Create new dyadic spatial proximity network
    sp_graph <- graph_from_adjacency_matrix(sp_net,mode="max")
    
  }
  
  #Create 1st higher-order layer (also referred to as the 'intermediate' layer (after Grueter et al. 2020 TREE))
  #Each maximal clique using the initial spatial proximity threshold is assigned as a hyperedge
  subgroups <- max_cliques(sp_graph)
  
  sp_hyp1 <- matrix(0, nrow = nrow(ind_data), ncol = length(subgroups))
  
  for(i in 1:length(subgroups)) {
    sp_hyp1[c(as.vector(subgroups[[i]])),i] <- 1 
  }
  
  sp_dyNet1 <- get_network(association_data = t(sp_hyp1), data_format = "GBI", association_index = "SRI")
  
  #Create higher-order layer 2 (also referred to as the 'upper' layer)
  thresh2 <- r[2]
  
  #Create matrix of subgroup membership
  sp_mat2 <- dist_mat < thresh2
  sp_net2 <- apply(sp_mat2, 2, as.numeric)
  sp_graph2 <- graph_from_adjacency_matrix(sp_net2, mode="max")
  
  #Create hypernetwork as above
  subgroups2 <- max_cliques(sp_graph2)
  sp_hyp2 <- matrix(0, nrow = nrow(ind_data), ncol = length(subgroups2))
  for(i in 1:length(subgroups2)) {
    sp_hyp2[c(as.vector(subgroups2[[i]])),i] <- 1 
  }
  
  sp_dyNet2 <- get_network(association_data = t(sp_hyp2), data_format = "GBI", association_index = "SRI")
  
  #Create higher-order layer 3 (also referred to as the 'apex' layer)
  thresh3 <- r[3]
  
  #Create matrix of subgroup membership
  sp_mat3 <- dist_mat < thresh3
  sp_net3 <- apply(sp_mat3, 2, as.numeric)
  sp_graph3 <- graph_from_adjacency_matrix(sp_net3, mode="max")
  
  subgroups3 <- max_cliques(sp_graph3)
  sp_hyp3 <- matrix(0, nrow = nrow(ind_data), ncol = length(subgroups3))
  
  for(i in 1:length(subgroups3)) {
    sp_hyp3[c(as.vector(subgroups3[[i]])),i] <- 1 
  }
  
  sp_dyNet3 <- get_network(association_data = t(sp_hyp3), data_format = "GBI", association_index = "SRI")
  
  return(list(sp_dyNet1, sp_hyp1, sp_dyNet2, sp_hyp2, sp_dyNet3, sp_hyp3))
}

hyperNetwork_diffusion <- function(ind_data, netList, network, informedNodes, domValues, resources, groupAdjustment = FALSE, 
                                   resources = resources) {
  
  produceTemp <- rep(0, nrow(ind_data))
  newLearner <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    focalHyperedge <- ifelse(length(which(netList[[network]][i,]>0)) > 1, sample(which(netList[[network]][i,]>0),1), which(netList[[network]][i,]>0))
    focalHyperedgeMembs <- which(netList[[network]][,focalHyperedge]>0)
    currentInformed <- focalHyperedgeMembs[focalHyperedgeMembs %in% informedNodes]
    
    if(length(currentInformed)==1){
      demons <- currentInformed
    } else{
      demons <- unique(sample(currentInformed, resources, replace = TRUE, 
                              prob = domValues[currentInformed]))
    }
    
    produceTemp[demons] <- 1
    currentUninformed <- focalHyperedgeMembs[!(focalHyperedgeMembs %in% currentInformed)]
    gA <- ifelse(groupAdjustment, length(focalHyperedgeMembs), 1)
    learningOutcomes <- runif(length(currentUninformed)) < (1 - (1 - social_trans/gA)^length(demons))
    
    if(length(currentUninformed[learningOutcomes]) > 0) {
      newLearner[currentUninformed[learningOutcomes]] <- 1
    }
  }
  
  return(list(produceTemp,newLearner))
}

dyadic_diffusion <- function(ind_data, netList, network, informedNodes, domValues, resources, groupAdjustment = FALSE, 
                             resources = resources) {
  
  produceTemp <- rep(0, nrow(ind_data))
  focalGraph <- graph_from_adjacency_matrix(netList[[network]], weighted = TRUE, mode = c("undirected"))
  
  for(i in informedNodes) {
    focalNeighbors <- as.vector(neighbors(focalGraph, v = i))
    currentInformed <- c(i, focalNeighbors[focalNeighbors %in% informedNodes])
    
    if(length(currentInformed) == 1) {
      demons <- currentInformed
    } else {
      demons <- unique(sample(currentInformed, resources, replace = TRUE, 
                              prob = domValues[currentInformed] * 
                                c(1,netList[[network]][i,focalNeighbors[focalNeighbors %in% informedNodes]])))
    }
    produceTemp[demons] <- 1
  }
  
  if(groupAdjustment) {
    acqProb <- sapply(1:nrow(ind_data), function(x) 1-prod(1-(produceTemp*netList[[network]][x,]*social_trans)/(rowSums(netList[[network]]))))
                                                               #1+sum(ifelse(ind_data$informed==1,0,1)*netList[[network]][x,])
  } else{
    acqProb <- sapply(1:nrow(ind_data), function(x) 1-prod(1-produceTemp*netList[[network]][x,]*social_trans))
  }
  newLearner <- ifelse(acqProb > runif(nrow(ind_data), min = 0, max = 1), 1, 0)
  newLearner[informedNodes] <- 0
  
  
  return(list(produceTemp,newLearner))
}