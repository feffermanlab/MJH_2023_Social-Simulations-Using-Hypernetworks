##Generate individuals in family/kinship groups and spatially embed them in x-y grid
generate_population <- function(n_families, meanFamilySize, clustering = NULL, nInitInformed, clusterByFamily) {
  #Calculate family size
  family_size <- rpois(n_families, meanFamilySize)
  
  #Calculate number of individuals
  n_indivs <- sum(family_size)
  
  #Assign IDs
  id <- seq(1, n_indivs, 1)
  
  #Assign Family IDs
  family <- c(rep(1:n_families, family_size))
  
  # #Assign dominance ranks
  # domRank <- runif(n_indivs, min = 0, max = 1)
  
  #Assign initial production probability for each individual
  initProduce <- produce <- runif(n_indivs, min = 0, max = 1)
  
  #Create individual dataframe
  ind_data<-data.frame(id, family, initProduce, produce)
  
  #Indicate for each individual the size of the kin group it belongs to
  ind_data$familySize <- sapply(1:nrow(ind_data), function(x) family_size[ind_data$family[x]])
  
  if(clusterByFamily == TRUE) {
  #Family coordinates
  family_x<-runif(n_families, 0, 1)
  family_y<-runif(n_families, 0, 1)
  
  #Assign individual coordinates
  ind_x <- ind_y <- rep(NA, nrow(ind_data))
  for(i in 1:nrow(ind_data)){
    ind_x[i] <- rnorm(1, family_x[ind_data$family[i]], clustering)
    ind_y[i] <- rnorm(1, family_y[ind_data$family[i]], clustering)
  }
  
  # #Rescale ind_x and ind_y to be between 0.0051 and 1
  # #0.0051 is selected as a minimum so that the round function will produce a value of at least 1 when a coordinate is multiplied by 100 (see resource matrix below)
  ind_data$ind_x <- ind_x
  ind_data$ind_y <- ind_y
  # ind_data$ind_x <- (ind_x - min(ind_x))/(max(ind_x) - min(ind_x)) * (1 - 0.0051) + 0.0051
  # ind_data$ind_y <- (ind_y - min(ind_y))/(max(ind_y) - min(ind_y)) * (1 - 0.0051) + 0.0051
  
  # #Assign colours to families for a plot
  # fam_cols <- sample(viridis::viridis(n_families), n_families, replace=FALSE)
  # 
  # #Can plot individual locations to check the algorithm is working OK
  # plot(ind_x, ind_y, col=fam_cols[ind_data$family], pch=16, cex=2)
  
  } else {
    if(clusterByFamily == FALSE) {
      #Assign individual coordinates
      ind_x <- ind_y <- rep(NA, nrow(ind_data))
      for(i in 1:nrow(ind_data)){
        ind_x[i] <- runif(1, 0, 1)
        ind_y[i] <- runif(1, 0, 1)
      }
      
      # #Rescale ind_x and ind_y to be between 0.0051 and 1
      # #0.0051 is selected as a minimum so that the round function will produce a value of at least 1 when a coordinate is multiplied by 100 (see resource matrix below)
      ind_data$ind_x <- ind_x
      ind_data$ind_y <- ind_y
    }
  }
  #Create an initial informed individual
  ind_data$informed <- 0
  ind_data$informed[sample(1:nrow(ind_data), nInitInformed, replace = FALSE)] <- 1
  
  #Identify initial informed nodes so that initial seed can be denoted
  informedNodes <- ind_data[which(ind_data$informed == 1),]$id
  ind_data$initInformed <- 0
  ind_data$initInformed[informedNodes] <- 1
  ind_data$acqTime <- 0
  
  return(ind_data)
}

##Create spatial proximity network and hypergraphs
generate_latent_space_multilayer_hypergraph_OLD <- function(ind_data, r) {
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
  
  return(list(sp_net, sp_hyp1, sp_hyp2, sp_hyp3))
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
  
  return(list(sp_net, sp_hyp1, sp_hyp2, sp_hyp3))
}

adjust_production_probability <- function(ind_data, netList, rule) {
  
  produceTemp <- ind_data$produce
  
  #Identify current set of informed nodes
  informedNodes <- ind_data[which(ind_data$informed == 1),]$id
  
  if(rule == "deltaAll") {
    #Informed individuals next adjust their current production probability to converge with that of other informed individuals in their intermediate-level hyperedges
    #delta all (general); see Sisk et al. (in prep)
    for(i in informedNodes) {
      adjustTemp <- mean(sapply(which(netList[[2]][i,] > 0), 
                                function(x) ifelse(sum(ind_data[which(ind_data$informed == 1),]$id %in% which(netList[[2]][,x] > 0)) <= 1, 0,
                                                   m * sum(sapply(which(netList[[2]][,x] > 0)[!(which(netList[[2]][,x] > 0) %in% i) &
                                                                                                (which(netList[[2]][,x] > 0) %in% informedNodes)],
                                                                  function(y) ind_data$produce[y] - ind_data$produce[i]))
                                )))
      produceTemp[i] <- min(1, max(0, ind_data$produce[i] + adjustTemp))
    }
  }
  
  if(rule == "deltaSmartest") {
    for(i in informedNodes) {
      #The following implements a version of the delta-smartest learning rule
      adjustTemp <- mean(sapply(which(netList[[2]][i,] > 0), 
                                function(x) ifelse(sum(ind_data[which(ind_data$informed == 1),]$id %in% which(netList[[2]][,x] > 0)) <= 1, 0,
                                                   m * (max(sapply(which(netList[[2]][,x] > 0)[which(netList[[2]][,x] > 0) %in% informedNodes],
                                                                   function(y) ind_data$produce[y])) - ind_data$produce[i])
                                )))
      produceTemp[i] <- min(1, max(0, ind_data$produce[i] + adjustTemp))
    }
  }
  
  if(rule == "dampMin") {
    for(i in informedNodes) {
      #The following is an adaptation of the dampening minimum rule
      adjustTemp <- sapply(which(netList[[2]][i,] > 0), 
                           function(x) ifelse(sum(ind_data[which(ind_data$informed == 1),]$id %in% which(netList[[2]][,x] > 0)) <= 1, ind_data$produce[i],
                                              (1/(sum(ind_data[which(ind_data$informed == 1),]$id %in% which(netList[[2]][,x] > 0)))) * (
                                                sum(sapply(which(netList[[2]][,x] > 0)[which(netList[[2]][,x] > 0) %in% informedNodes],
                                                           function(y) ind_data$produce[y])) - 
                                                  (m * (max(sapply(which(netList[[2]][,x] > 0)[which(netList[[2]][,x] > 0) %in% informedNodes],
                                                                   function(y) ind_data$produce[y])) - min(sapply(which(netList[[2]][,x] > 0)[which(netList[[2]][,x] > 0) %in% informedNodes],
                                                                                                                  function(y) ind_data$produce[y])))))
                           ))
      adjustTemp <- mean(ifelse(adjustTemp > ind_data$produce[i], adjustTemp, ind_data$produce[i]))
      produceTemp[i] <- min(1, max(0, adjustTemp))
    }
  }
  
  if(rule == "deltaCombined") {
    for(i in informedNodes) {
      adjustTemp <- mean(sapply(which(netList[[2]][i,] > 0), function(x) 
        ifelse(sum(ind_data[which(ind_data$informed == 1),]$id %in% which(netList[[2]][,x] > 0)) <= 1, 0,
               ifelse(sum(ind_data[which(ind_data$informed == 1),]$id %in% which(netList[[2]][,x] > 0)) / length(which(netList[[2]][,x] > 0)) < 0.5,
                      m * (min(sapply(which(netList[[2]][,x] > 0)[which(netList[[2]][,x] > 0) %in% informedNodes],
                                      function(y) ind_data$produce[y])) - ind_data$produce[i]),
                      m * (max(sapply(which(netList[[2]][,x] > 0)[which(netList[[2]][,x] > 0) %in% informedNodes],
                                      function(y) ind_data$produce[y])) - ind_data$produce[i]))
        )))
      produceTemp[i] <- min(1, max(0, ind_data$produce[i] + adjustTemp))
    }
  }
  return(produceTemp)
}

produceByDomRank <- function(ind_data, netList, domSteepness, netType) {
  
  produceTemp <- ind_data$produce
  
  #Identify current set of informed nodes
  informedNodes <- ind_data[which(ind_data$informed == 1),]$id  
  
  if(netType == "higherOrder") {
  for(i in informedNodes) {
    adjustTemp <- mean(
      sapply(which(netList[[2]][i,] > 0), function(x) sum(sapply(which(netList[[2]][,x] > 0)[which(!(which(netList[[2]][,x]>0) %in% i))], 
                                                                                  function(y) ifelse(domSteepness * (ind_data$domRank[i]/(ind_data$domRank[i] + ind_data$domRank[y])) > runif(1,0,1), 1, 0))) 
                              / (length(which(netList[[2]][,x]>0)) - 1))
      )
    
    produceTemp[i] <- adjustTemp
  }
  } else{
      if(netType == "dyadic") {
        for(i in informedNodes) {
        adjustTemp <- sum(sapply(which(netList[[1]][i,]> 0), function(y) ifelse(domSteepness * (ind_data$domRank[i]/(ind_data$domRank[i] + ind_data$domRank[y])) > runif(1,0,1), 1, 0)))/length(which(netList[[1]][i,]>0))
        produceTemp[i] <- adjustTemp
        }
      }
    }
  return(produceTemp)
}

domResponse_Linear_HE <- function(ind_data, netList, radius, informedNodes, domValues) {
  
  produceTemp <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    adjustTemp <- mean(
      sapply(which(netList[[radius]][i,] > 0), function(x) 
        1 - (sum(sapply(which(netList[[radius]][,x] > 0)[which(!(which(netList[[radius]][,x]>0) %in% i))], 
                   function(y) ifelse(domValues[i]/(domValues[i] + domValues[y]) > runif(1,0,1), 0, 1))) /
          (length(which(netList[[radius]][,x]>0)) - 1))
        )
    )
    
    produceTemp[i] <- adjustTemp
  }
  
  return(produceTemp)
}

domResponse_Sigmoid1_HE <- function(ind_data, netList, radius, informedNodes, domValues) {
  
  produceTemp <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    adjustTemp <- mean(
      sapply(which(netList[[radius]][i,] > 0), function(x) 
        1 / (1 + exp(sum(sapply(which(netList[[radius]][,x] > 0)[which(!(which(netList[[radius]][,x]>0) %in% i))], 
                        function(y) ifelse(domValues[i]/(domValues[i] + domValues[y]) > runif(1,0,1), 0, 1))) - 
                       ((length(which(netList[[radius]][,x]>0)) - 1)/2))
               )
      )
    )
    
    produceTemp[i] <- adjustTemp
  }
  
  return(produceTemp)
}


domResponse_Sigmoid2_HE <- function(ind_data, netList, radius, informedNodes, domValues) {
  
  produceTemp <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    adjustTemp <- mean(
      sapply(which(netList[[radius]][i,] > 0), function(x) 
        1 / (1 + exp(0.75 * sum(sapply(which(netList[[radius]][,x] > 0)[which(!(which(netList[[radius]][,x]>0) %in% i))], 
                                function(y) ifelse(domValues[i]/(domValues[i] + domValues[y]) > runif(1,0,1), 0, 1))) - 
                       ((length(which(netList[[radius]][,x]>0)) - 1)/5))
        )
      )
    )
    
    produceTemp[i] <- adjustTemp
  }
  
  return(produceTemp)
}


domResponse_Linear_dyadic <- function(ind_data, netList, radius, informedNodes, domValues) {
  
  produceTemp <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    adjustTemp <- 1 - (sum(sapply(which(netList[[radius]][i,]> 0), function(y) ifelse(domValues[i]/(domValues[i] + domValues[y]) > runif(1,0,1), 0, 1)))/length(which(netList[[radius]][i,]>0)))
    produceTemp[i] <- adjustTemp
  }
  
  return(produceTemp)
}

domResponse_Sigmoid1_dyadic <- function(ind_data, netList, radius, informedNodes, domValues) {
  
  produceTemp <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    adjustTemp <- 1 / (1 + exp(sum(sapply(which(netList[[radius]][i,]> 0), function(y) 
      ifelse(domValues[i]/(domValues[i] + domValues[y]) > runif(1,0,1), 0, 1))) - 
        (length(which(netList[[radius]][i,]>0))/2)))
    produceTemp[i] <- adjustTemp
  }
  
  return(produceTemp)
}


domResponse_Sigmoid2_dyadic <- function(ind_data, netList, radius, informedNodes, domValues) {
  
  produceTemp <- rep(0, nrow(ind_data))
  
  for(i in informedNodes) {
    adjustTemp <- 1 / (1 + exp(0.75 * sum(sapply(which(netList[[radius]][i,]> 0), function(y) 
      ifelse(domValues[i]/(domValues[i] + domValues[y]) > runif(1,0,1), 0, 1))) - 
        (length(which(netList[[radius]][i,]>0))/5)))
    produceTemp[i] <- adjustTemp
  }
  
  return(produceTemp)
}