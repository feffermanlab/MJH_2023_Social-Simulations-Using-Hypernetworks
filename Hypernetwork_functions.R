
#Create a hypernetwork similar to an Erdos-Renyi graph
#Each vertex, v, has a probability, p, of being in each of m hyperedges
#See Aksoy et al., 2020, EPJ Data Science, 9: 16 for more details
generate_ER_hypergraph <- function(v, m, p) {
  hyperEdgeList <- vector("list", m)
  for(i in 1:v) {
    for(j in 1:m) {
      if(runif(1, min = 0, max = 1) <= p) {
        hyperEdgeList[[j]] <- c(hyperEdgeList[[j]], i)
      }
    }
  }
  
  if(sum(sapply(hyperEdgeList, is.null)) > 0) {
    hyperEdgeList <- hyperEdgeList[-which(sapply(hyperEdgeList, is.null))]
  }
  
  for(i in 1:v) {
    if(sum(sapply(hyperEdgeList, function(x) i %in% x)) == 0) {
      hyperEdgeList[length(hyperEdgeList) + 1] <- i
    }
  }
  
  hyperEdgeList <- hyperEdgeList[order(sapply(hyperEdgeList, '[[', 1))]
  hyperEdgeList <- hyperEdgeList[order(sapply(hyperEdgeList, length))]
  
  return(hyperEdgeList)
}

create_population_data <- function(N, initialAges, maxIndivs = 10000, ageBias, selectGradient) {
  dataTemp <- data.frame("ID" = seq(1:maxIndivs), 
                         "Age" = c(initialAges, rep(0, maxIndivs - N)),
                         "Parent" = 0,
                         "GSPref" = c(sample(seq(from = 2, to = 10), size = N, replace = TRUE), rep(0, maxIndivs - N)),
                         "AgeBias" = ageBias,
                         "selectGradient" = selectGradient,
                         "Alive" = c(rep("Y", N), rep("NB", maxIndivs - N)))
}

eidx <- function(t) {
  paste("e", t, sep="")
}

get_dual_hypergraph <- function(hypNet, popData) {
  hyperEdgeList <- hypNet
  dualHyperEdgeList <- vector("list", nrow(popData[which(popData$Alive == "Y"),]))
  names(dualHyperEdgeList) <- c(popData$ID[which(popData$Alive == "Y")])
  for(i in popData$ID[which(popData$Alive == "Y")]) {
    edgeTemp <- c(which(sapply(hyperEdgeList, function(x) i %in% x)))
    dualHyperEdgeList[[as.character(i)]] <- c(edgeTemp)
  }
  return(dualHyperEdgeList)
}

get_s_line_graph <- function(hypNet, size, vertexNames, mode){
  if(mode == "incidence") {
    matrixTemp <- matrix(0, nrow = ncol(hypNet), ncol = ncol(hypNet))
    rownames(matrixTemp) <- vertexNames
    colnames(matrixTemp) <- vertexNames
    hyperEdgeSizes <- sapply(seq(from = 1, to = ncol(hypNet)), function(x) sum(hypNet[,x] > 0))
    hyperEdgeIndices <- which(hyperEdgeSizes >= size)
    #What is the purpose of the following line? Is it supposed to be sizes?
    #hyperEdgeIndices <- hyperEdgeIndices[which(hyperEdgeIndices != ncol(hypNet))]
    if(length(hyperEdgeIndices) > 1) {
    for(i in hyperEdgeIndices) {
      refVect <- as.vector(which(hypNet[,i]>0))
      focalCol <- hypNet[,i]
      indices <- hyperEdgeIndices[which(hyperEdgeIndices != i)]
      matrixTemp[i,c(indices)] <- sapply(indices, 
                         function(x)
                           ifelse(
                           length(focalCol[c(intersect(refVect, 
                                                      as.vector(which(hypNet[,x]>0))))]) >= size, 1, 0))
    }
    matrixTemp[lower.tri(matrixTemp)] <- t(matrixTemp)[lower.tri(matrixTemp)]
    }
  }
  else{
    if(mode == "edgeList") {
  hyperEdgeList <- hypNet
  matrixTemp <- matrix(0, nrow = length(hyperEdgeList), ncol = length(hyperEdgeList))
  rownames(matrixTemp) <- names(hyperEdgeList)
  colnames(matrixTemp) <- names(hyperEdgeList)
  
  for(i in 1:length(hyperEdgeList)) {
    for(j in 1:length(hyperEdgeList)) {
      if(i != j){
        if(length(intersect(hyperEdgeList[[i]], hyperEdgeList[[j]])) >= size) {
          matrixTemp[i,j] <- 1
        }
      }
    }
  }
  matrixTemp <- matrixTemp[c(as.character(sort(as.integer(rownames(matrixTemp)))))
                           ,c(as.character(sort(as.integer(colnames(matrixTemp)))))]
    }
  }
  return(matrixTemp)
}

get_weighted_graph <- function(hypNet){
  hyperEdgeList <- hypNet
  matrixTemp <- matrix(0, nrow = length(hyperEdgeList), ncol = length(hyperEdgeList))
  rownames(matrixTemp) <- names(hyperEdgeList)
  colnames(matrixTemp) <- names(hyperEdgeList)
  
  for(i in 1:length(hyperEdgeList)) {
    for(j in 1:length(hyperEdgeList)) {
      if(i != j){
        if(length(intersect(hyperEdgeList[[i]], hyperEdgeList[[j]])) > 0) {
          matrixTemp[i,j] <- length(intersect(hyperEdgeList[[i]], hyperEdgeList[[j]]))
        }
      }
    }
  }
  matrixTemp <- matrixTemp[c(as.character(sort(as.integer(rownames(matrixTemp)))))
                           ,c(as.character(sort(as.integer(colnames(matrixTemp)))))]
  return(matrixTemp)
}

select_partners <- function(prefMatrix, popData, t) {
  livePop <- popData[which(popData$Alive == "Y"),]
  groupList <- vector("list", nrow(livePop))
  pairList <- matrix(0, nrow = nrow(livePop), ncol = 2)
  selectOrder <- sample(livePop$ID, nrow(livePop), 
                        prob = livePop$Age, 
                        replace = FALSE)
  prefMatrixTemp <- prefMatrix[[t]]
  for(i in as.vector(which(rowSums(prefMatrixTemp) == 0))) {
    prefMatrixTemp[i,] <- 1
    prefMatrixTemp[i, i] <- 0
  }
  for(i in selectOrder) {
    if(length(groupList[sapply(groupList, function(x) i %in% x)]) == 0) {
      currentAge <- livePop$Age[which(livePop$ID == i)]
      selectivity <- 1 + 1/(1 + exp(-(livePop$selectGradient[which(livePop$ID == i)]) * currentAge + 5))
      associates <- as.integer(colnames(prefMatrixTemp)[prefMatrixTemp[as.character(i),] > 0])
      associateData <- popData[which(popData$ID %in% associates),]
      associateData$PrevInt <- sapply(associates, function(x) prefMatrixTemp[as.character(i), as.character(x)])
      associateData$PrevIntSelective <- associateData$PrevInt ^ selectivity
      associateData$RelativePrevInt <- associateData$PrevIntSelective/sum(associateData$PrevIntSelective)
      associateData$RelativeAge <- associateData$Age/sum(associateData$Age)
      associateData$AgeDifference <- abs(currentAge - associateData$Age) + 1
      ifelse(nrow(associateData) == 1, 
             associateData$AgeDiffRev <- 1, 
      associateData$AgeDiffRev <- 1 - (associateData$AgeDifference/sum(associateData$AgeDifference)))
      associateData$RelativeAgeDiff <- associateData$AgeDiffRev/sum(associateData$AgeDiffRev)
      associateData$PsUnadjust <- ((livePop$AgeBias[which(livePop$ID == i)] * associateData$RelativeAge) + 
                                     ((1 - livePop$AgeBias[which(livePop$ID == i)]) * associateData$RelativeAgeDiff))
      associateData$PsAdjust <- associateData$PsUnadjust/sum(associateData$PsUnadjust)
      associateData$Ps <- associateData$PsAdjust * associateData$RelativePrevInt
      preferenceOrder <- as.integer(sample(as.character(associateData$ID), nrow(associateData), prob = associateData$Ps, replace = FALSE))
      
      repeat{
        j = preferenceOrder[1]
        if(length(groupList[sapply(groupList, function(x) j %in% x)]) == 0) {
          groupList[[1 + (length(groupList) - length(groupList[which(sapply(groupList, is.null))]))]] <- c(i, j)
          pairList[which(livePop$ID == i),] <- c(i,j)
          break
        } else{
          if(length(groupList[sapply(groupList, function(x) j %in% x)][[1]]) <
             min(
               c(popData$GSPref[which(popData$ID == i)], 
                 popData$GSPref[which(popData$ID %in% groupList[sapply(groupList, function(x) j %in% x)][[1]])]))) {
            groupList[[which(sapply(groupList, function(x) j %in% x))[1]]] <- c(groupList[[which(sapply(groupList, function(x) j %in% x))[1]]], i)
            pairList[which(livePop$ID == i),] <- c(i,j)
            break
          }
        }
        preferenceOrder <- preferenceOrder[-1]
        if(length(preferenceOrder) == 0){
          groupList[[1 + (length(groupList) - length(groupList[which(sapply(groupList, is.null))]))]] <- c(i)
          break
        }
      }
    }
  }
  groupList <- groupList[-which(sapply(groupList, is.null))]
  pairList <- pairList[pairList[,1] > 0, ]
  return(list(groupList, pairList))
}

generate_age_structure <- function(N, maxIter) {
  dataTemp <- data.frame("ID" = seq(1:N),
                         "Age" = 1,
                         "Alive" = "Y")
  iteration = 1
  repeat{
    if(iteration ==  maxIter){
      break
    }
    
    dataTemp$Alive <- sapply(dataTemp$Age, function(x) sample(c("Y","N"), 1, prob = c(
      1 - ((10^-1)*(1/(1+exp((-0.1*x)+5)))), 
      ((10^-1)*(1/(1+exp((-0.1*x)+5)))))))
    
    if(sum(dataTemp$Alive == "N") > 0){
      for(i in 1:sum(dataTemp$Alive == "N")) { 
        maxID = max(dataTemp$ID)
        newIndData <- data.frame("ID" = maxID + 1, 
                                 "Age" = 0,
                                 "Alive" = "Y")
        dataTemp <- rbind(dataTemp, newIndData)
      }
    }
    
    dataTemp <- dataTemp[which(dataTemp$Alive == "Y"),]
    dataTemp$Age <- dataTemp$Age + 1
    iteration <- iteration + 1
  }
  return(dataTemp$Age)
}

# generate_full_knowledge_set <- function(numDomains, numLevels, complexity) {
#   domains <- letters[1:numDomains]
#   
#   #May want to modify this eventually to allow for different number of levels in each domain
#   levels <- seq(from = 1, to = numLevels)
#   
#   simpleKnowledgeSet <- list()
#   
#   for(i in domains){
#     for(j in levels){
#       simpleKnowledgeSet <- c(simpleKnowledgeSet, 
#                               paste(i, j, sep = ""))
#     }
#   }
#   knowledgeSet <- simpleKnowledgeSet
#   
#   if(complexity > 1) {
#     complexKnowledgeSet <- list()
#     for(c in 2:complexity){
#       for(i in domains){
#         kTemp <- grep(i, simpleKnowledgeSet, value = TRUE)
#         complexKnowledgeSet <- c(complexKnowledgeSet, c(combn(kTemp, c, simplify = FALSE, FUN = paste, collapse = '')))
#       }
#     }
#     knowledgeSet <- c(simpleKnowledgeSet, complexKnowledgeSet)
#   }
#   
#   return(knowledgeSet)
# }
# 
# generate_initial_pop_knowledge <- function(popData, knowledgeSet, asoc, simpleKnowledge = TRUE){
#   popKnowledge <- vector("list", nrow(popData))
#   
#   #Works for simple knowledge; will need to decide how to handle complex knowledge (e.g., "a1 a2")
#   if(simpleKnowledge == TRUE) {
#     simpleKnowl <- knowledgeSet[sapply(knowledgeSet, function(x) length(x) == 1)]
#     popAges <- popData$Age[which(popData$Alive == "Y")]
#     numDiscoveries <- rep(0, length(popAges))
#     for(i in 1:length(popAges)) {
#       numDiscoveries[i] <- sum(sample(c(0,1), popAges[i], replace = TRUE, prob = c(1 - asoc, asoc)))
#     }
#     for(i in 1:length(numDiscoveries)) {
#       popKnowledge[i] <- list(sample(simpleKnowl, numDiscoveries[i], replace = FALSE))
#     }
#   } else{
#     fullKnowl <- knowledgeSet
#     popAges <- popData$Age[which(popData$Alive == "Y")]
#     numDiscoveries <- rep(0, length(popAges))
#     for(i in 1:length(popAges)) {
#       numDiscoveries[i] <- sum(sample(c(0,1), popAges[i], replace = TRUE, prob = c(1 - asoc, asoc)))
#     }
#     for(i in 1:length(numDiscoveries)) {
#       popKnowledge[i] <- list(sample(fullKnowl, numDiscoveries[i], replace = FALSE))
#     }
#   }
#   return(popKnowledge)
# }

initialize_preference_matrix <- function(N, maxT, currentPartners) {
  matrixList <- lapply(1:maxT, matrix, data = 0, nrow = N, ncol = N)
  currentMatrix <- matrixList[[1]]
  rownames(currentMatrix) <- 1:N
  colnames(currentMatrix) <- 1:N
  for(i in 1:length(currentPartners)) {
    if(length(currentPartners[[i]]) > 1) {
      combnMatrix <- combn(currentPartners[[i]], 2)
      for(j in 1:dim(combnMatrix)[2]) {
        currentMatrix[as.character(combnMatrix[,j][1]), as.character(combnMatrix[,j][2])] <- 
          1 + currentMatrix[as.character(combnMatrix[,j][1]), as.character(combnMatrix[,j][2])]
        currentMatrix[as.character(combnMatrix[,j][2]), as.character(combnMatrix[,j][1])] <- 
          1 + currentMatrix[as.character(combnMatrix[,j][2]), as.character(combnMatrix[,j][1])]
      }
    }
  }
  matrixList[[1]] <- currentMatrix
  return(matrixList)
}

# get_hyperedge_knowledge <- function(hyperEdges, popKnowl) {
#   hypEKnowledge <- data.frame("edgeID" = 1:length(hyperEdges), 
#                               "KTotal" = 0,
#                               "KDivers" = 0)
#   for(i in 1:length(hyperEdges)) {
#     kTemp <- unlist(sapply(hyperEdges[[i]], function(x) unlist(popKnowl[[x]])))
#     hypEKnowledge$KTotal[hypEKnowledge$edgeID == i] <- length(kTemp)
#     hypEKnowledge$KDivers[hypEKnowledge$edgeID == i] <- length(unique(kTemp))
#   }
#   return(hypEKnowledge)
# }

# learn_new_knowledge <- function(hyperEdgeKnowledge, popData, popKnowl, hyperEdges, selectMethod, kThreshold) {
#   livePop <- popData[which(popData$Alive == "Y"),]
#   newKnowledge <- vector("list", nrow(livePop))
#   names(newKnowledge) <- c(livePop$ID)
#   if(selectMethod == "total"){
#     maxTotal <- rep(0, length(livePop$ID))
#   for(i in 1:length(livePop$ID)) {
#     currentEdges <- c(which(sapply(hyperEdges, function (x) livePop$ID[i] %in% x)))
#       maxTotal[i] <- max(hyperEdgeKnowledge$KTotal[which(hyperEdgeKnowledge$edgeID %in% currentEdges)])
#   }
#     for(i in which(maxTotal > 0)){
#       #if(maxTotal > 0) {
#         currentEdges <- c(which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x)))
#         selectedEdge <- hyperEdgeKnowledge$edgeID[which(hyperEdgeKnowledge$edgeID %in% currentEdges & 
#                                                           hyperEdgeKnowledge$KTotal == maxTotal[i])]
#         #if(length(selectedEdge) > 1) {
#           selectedEdge <- as.integer(sample(as.character(selectedEdge), size = 1))
#         #}
#         hyperEdgeKnowlVariants <- unlist(sapply(hyperEdges[[selectedEdge]], function(x) unlist(popKnowl[[x]])))
#         variantData <- data.frame("variantID" = sort(unique(hyperEdgeKnowlVariants)), 
#                                   "propKnowledgeable" = 0)
#         for(r in 1:nrow(variantData)) {
#           variantData[r,"propKnowledgeable"] <- sum(hyperEdgeKnowlVariants == variantData[r, "variantID"])/length(hyperEdges[[selectedEdge]])
#         }
#         learnableVariants <- variantData[which(variantData$propKnowledgeable > kThreshold),]
#         if(nrow(learnableVariants) > 0) {
#           learnedVariants <- learnableVariants$variantID[which(!(learnableVariants$variantID %in% unlist(popKnowl[[livePop$ID[i]]])))]
#           if(length(learnedVariants) > 0){
#             newKnowledge[as.character(livePop$ID[i])] <- list(learnedVariants)
#           }
#           }
#       }
#     } 
#     else{
#       if(selectMethod == "diversity") {
#         uniqueTotal <- max(hyperEdgeKnowledge$KDivers[which(hyperEdgeKnowledge$edgeID %in% currentEdges)])
#         selectedEdge <- hyperEdgeKnowledge$edgeID[which(hyperEdgeKnowledge$edgeID %in% currentEdges & 
#                                                           hyperEdgeKnowledge$KDivers == uniqueTotal)]
#         if(length(selectedEdge) > 1) {
#           selectedEdge <- sample(selectedEdge, n = 1)
#         } 
#       }
#     }
#   return(newKnowledge)
# }

# update_knowledge <- function(popKnowl, newKnowl) {
#   newKnowl[sapply(newKnowl, is.null)] <- NULL
#   for(r in 1:length(newKnowl)) {
#     popKnowl[[as.integer(names(newKnowl)[[r]])]] <- append(popKnowl[[as.integer(names(newKnowl)[[r]])]], as.list(newKnowl[[r]]))
#   }
#   return(popKnowl)
# }

identify_demographic_changes <- function(popData) {
  livePop <- popData[which(popData$Alive == "Y"),]
  livePop$Alive <- sapply(livePop$Age, function(x) sample(c("Y","N"), 1, prob = c(
    1 - (0.1/(1+exp((-0.1*x)+5))), 
    0.1/(1+exp((-0.1*x)+5))
  )
  )
  )
  numDeaths <- sum(livePop$Alive == "N")
  deathIdentities <- livePop$ID[which(livePop$Alive == "N")]
  parentIdentities <- sample(livePop$ID, size = numDeaths, replace = FALSE, prob = ((10^-1)*(1/(1+exp((-0.1*livePop$Age)+5)))))
  return(list(deathIdentities, parentIdentities))
}

update_preference_matrix <- function(prefMatrices, popData, timeStep, demoChanges, currentPartners, inheritance) {
  currentMatrix <- prefMatrices[[timeStep]]
  for(i in 1:length(currentPartners)) {
    if(length(currentPartners[[i]]) > 1) {
    combnMatrix <- combn(currentPartners[[i]], 2)
    for(j in 1:dim(combnMatrix)[2]) {
      currentMatrix[as.character(combnMatrix[,j][1]), as.character(combnMatrix[,j][2])] <- 
        1 + currentMatrix[as.character(combnMatrix[,j][1]), as.character(combnMatrix[,j][2])]
      currentMatrix[as.character(combnMatrix[,j][2]), as.character(combnMatrix[,j][1])] <- 
        1 + currentMatrix[as.character(combnMatrix[,j][2]), as.character(combnMatrix[,j][1])]
    }
    }
  }
  if(length(demoChanges[[2]]) > 0) {
    maxID <- max(popData$ID[which(popData$Alive == "Y")])
    newIDs <- rep(0, length(demoChanges[[2]]))
    for(i in 1:length(demoChanges[[2]])) {
      newIDs[i] <- maxID + i
    }
    oldNames <- colnames(currentMatrix)
    newIDColumns <- matrix(0, nrow = nrow(currentMatrix), ncol = length(newIDs))
    colnames(newIDColumns) <- as.character(newIDs)
    rownames(newIDColumns) <- row.names(currentMatrix)
    currentMatrix <- cbind(currentMatrix, newIDColumns)
    newIDRows <- matrix(0, nrow = length(newIDs), ncol = ncol(currentMatrix))
    rownames(newIDRows) <- as.character(newIDs)
    colnames(newIDRows) <- colnames(currentMatrix)
    currentMatrix <- rbind(currentMatrix, newIDRows)
    
    if(inheritance == "random") {
      for(i in 1:length(newIDs)) {
        parent <- demoChanges[[2]][[i]]
        currentEdges <- c(which(sapply(currentPartners, function (x) parent %in% x)))
        numPartners <- length(unique(unlist(currentPartners[currentEdges])[which(unlist(currentPartners[currentEdges]) != parent)])) + 1
        possiblePartners <- as.integer(colnames(currentMatrix))[which(colnames(currentMatrix) != newIDs[i])]
        partners <- sample(possiblePartners, numPartners, replace = FALSE)
        for(j in partners) {
          currentMatrix[as.character(newIDs[i]), as.character(j)] <- 1
          currentMatrix[as.character(j), as.character(newIDs[i])] <- 1
        }
        currentMatrix[as.character(newIDs[i]), as.character(parent)] <- 1
        currentMatrix[as.character(parent), as.character(newIDs[i])] <- 1
      }
    } else{
      if(inheritance == "parental") {
        for(i in 1:length(newIDs)) {
          parent <- demoChanges[[2]][[i]]
          currentEdges <- c(which(sapply(currentPartners, function (x) parent %in% x)))
          partners <- unique(unlist(currentPartners[currentEdges]))
                             #[which(unlist(currentPartners[currentEdges]) != parent)])
          for(j in partners) {
            currentMatrix[as.character(newIDs[i]), as.character(j)] <- 1
            currentMatrix[as.character(j), as.character(newIDs[i])] <- 1
          }
          currentMatrix[as.character(newIDs[i]), as.character(parent)] <- 1
          currentMatrix[as.character(parent), as.character(newIDs[i])] <- 1
        } 
      }
    }
    
    for(i in demoChanges[[1]]) {
      currentMatrix <- currentMatrix[!rownames(currentMatrix) %in% as.character(i), !colnames(currentMatrix) %in% as.character(i)]
    }
  }
  return(currentMatrix) 
}


update_population_data <- function(popData, demoChanges, ageBias, selectGradient) {
  
  if(length(demoChanges[[2]]) > 0) {
    maxID <- max(popData$ID[which(popData$Alive == "Y")])
    newIDs <- rep(0, length(demoChanges[[2]]))
    for(i in 1:length(demoChanges[[2]])) {
      newIDs[i] <- maxID + i
    }
    newIndivs <- data.frame("ID" = newIDs,
                            "Age" = 1,
                            "Parent" = demoChanges[[2]],
                            "GSPref" = sample(seq(from = 2, to = 10), size = length(newIDs), replace = TRUE),
                            "AgeBias" = ageBias,
                            "selectGradient" = selectGradient,
                            "Alive" = "Y")
    for(i in demoChanges[[1]]) {
      popData$Alive[popData$ID == i] <- "N"
    }
    popData$Age[popData$Alive == "Y"] <- popData$Age[popData$Alive == "Y"] + 1
    popData[newIDs,] <- newIndivs
  } else{
    popData$Age[popData$Alive == "Y"] <- popData$Age[popData$Alive == "Y"] + 1
  }
  return(popData)
}

# generate_knowledge_for_new_indivs <- function(popData, knowledgeSet, asoc, simpleKnowledge = TRUE){
#   knowledgeTemp <- vector("list", length(popData$ID[which(popData$Age == 1)]))
#   
#   #Works for simple knowledge; will need to decide how to handle complex knowledge (e.g., "a1 a2")
#   if(simpleKnowledge == TRUE) {
#     simpleKnowl <- knowledgeSet[sapply(knowledgeSet, function(x) length(x) == 1)]
#     numDiscoveries <- sample(c(0,1), length(which(popData$Age == 1)), replace = TRUE, prob = c(1 - asoc, asoc))
#     for(i in 1:length(numDiscoveries)) {
#       knowledgeTemp[i] <- list(sample(simpleKnowl, numDiscoveries[i], replace = FALSE))
#     }
#   }
#   return(knowledgeTemp)
# }

# asocial_learning <- function(popKnowledge, popData, knowledgeSet, simpleKnowledge, asoc) {
#   if(simpleKnowledge == TRUE) {
#     numDiscoveries <- sample(c(0,1), length(which(popData$Age > 1 & popData$Alive == "Y")), replace = TRUE, prob = c(1-asoc, asoc))
#     livePopIDs <- popData$ID[popData$Alive == "Y" & popData$Age > 1]
#     for(i in which(numDiscoveries == 1)) {
#       knowledgePool <- which(!(knowledgeSet[sapply(knowledgeSet, function(x) length(x) == 1)] %in% unlist(popKnowledge[[livePopIDs[i]]])))
#       if(length(knowledgePool) > 0) {
#         newKnowl <- sample(knowledgePool, size = 1)
#         popKnowledge[[livePopIDs[i]]] <- append(popKnowledge[[livePopIDs[i]]], knowledgeSet[[newKnowl]])
#       }
#     }
#   }
#   return(popKnowledge)
# }

get_s_betweenness <- function(hypergraph, smax, vertexNames, mode) {
  matrixTemp <- matrix(0, nrow = ncol(hypergraph), ncol = smax)
  matrixTemp <- sapply(seq_along(1:smax), function(x) 
    matrixTemp[,x] <- as.vector(
      betweenness(
        graph_from_adjacency_matrix(
          get_s_line_graph(
            hypergraph, size = x, vertexNames = vertexNames, mode = mode), 
          mode = "undirected"), 
        directed = FALSE, normalized = TRUE)))
  rankedMatrix <- matrix(0, nrow = ncol(hypergraph), ncol= smax)
  rankedMatrix <- sapply(seq_along(1:smax), function(x) rankedMatrix[,x] <- 
                           rank(matrixTemp[,x]))
  integratedScores <- rank(-rowSums(rankedMatrix)/smax, ties.method = "average")
  return(list(matrixTemp, integratedScores))
}

#Does not adequately function right now due to s-line graphs becoming disconnected for higher values of s
get_s_eigen_centrality <- function(hypergraph, smax, vertexNames, mode) {
  matrixTemp <- matrix(0, nrow = ncol(hypergraph), ncol = smax)
  matrixTemp <- sapply(seq_along(1:smax), function(x) 
    matrixTemp[,x] <- as.vector(
      eigen_centrality(
        graph_from_adjacency_matrix(
          get_s_line_graph(
            hypergraph, size = x, vertexNames = vertexNames, mode = mode), 
          mode = "undirected", diag = FALSE), 
        directed = FALSE, scale = TRUE)$vector))
  rankedMatrix <- matrix(0, nrow = ncol(hypergraph), ncol= smax)
  rankedMatrix <- sapply(seq_along(1:smax), function(x) rankedMatrix[,x] <- 
                           rank(matrixTemp[,x]))
  integratedScores <- rank(-rowSums(rankedMatrix)/smax, ties.method = "average")
  return(list(matrixTemp, integratedScores))
}

get_s_harmonic_centrality <- function(hypergraph, smax, vertexNames, mode) {
  matrixTemp <- matrix(0, nrow = length(hypergraph), ncol = smax)
  matrixTemp <- sapply(seq_along(1:smax), function(x) 
    matrixTemp[,x] <- as.vector(
      harmonic_centrality(
        graph_from_adjacency_matrix(
          get_s_line_graph(
            hypergraph, size = x, vertexNames = vertexNames, mode = mode), 
          mode = "undirected"), 
        normalized = TRUE)))
  rankedMatrix <- matrix(0, nrow = length(hypergraph), ncol= smax)
  rankedMatrix <- sapply(seq_along(1:smax), function(x) rankedMatrix[,x] <- 
                           rank(matrixTemp[,x]))
  integratedScores <- rank(-rowSums(rankedMatrix)/smax, ties.method = "average")
  return(list(matrixTemp, integratedScores))
}

# get_number_of_variants_known  <- function(popData, popKnowledge, complexity = NULL) {
#   liveIDs <- popData$ID[which(popData$Alive == "Y")]
#   variantData <- rep(0, length(liveIDs))
#   ifelse(is.null(complexity), 
#          variantData <- sapply(liveIDs, 
#                                function(x) length(unlist(popKnowledge[[x]]))), 
#          variantData <- sapply(liveIDs, 
#                                function(x) ifelse(length(popKnowledge[[x]]) == 0, 0, 
#                                                   sum(sapply(popKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == complexity))
#                                )
#          )
#   )
#   return(variantData)
# }

# record_individual_summary_data <- function(popData, popKnowledge = NULL, knowledgeSet = NULL, smax, mode) {
#   if(is.null(knowledgeSet)) {    
#     data.temp <- data.frame(
#     "ID" = popData$ID[which(popData$Alive == "Y")], 
#     "timeStep" = t, 
#     "Age" = popData$Age[which(popData$Alive == "Y")],
#     "Parent" = popData$Parent[which(popData$Alive == "Y")],
#     "GSPref" = popData$GSPref[which(popData$Alive == "Y")],
#     "s1BC" = 0, "s2BC" = 0, "s3BC" = 0,
#     "siBC" = 0,
#     "s1D" = 0, "s2D" = 0, "s3D" = 0,
#     "siD" = 0, 
#     "s1HC" = 0, "s2HC" = 0, "s3HC" = 0,
#     "siHC" = 0, 
#     "degree" = 0,
#     "strength" = 0, 
#     "betweenness" = 0, 
#     "harmonicCent" = 0)
#   dualHyperNetwork <- get_dual_hypergraph(hypNet = currentPartners, popData = popData)
#   sBCScores <- get_s_betweenness(hypergraph = dualHyperNetwork, smax = smax, mode = mode)
#   sDegreeScores <- get_s_degree(hypergraph = dualHyperNetwork, smax = smax, mode = mode)
#   sHCScores <- get_s_harmonic_centrality(hypergraph = dualHyperNetwork, smax = smax, mode = mode)
#   data.temp[,c('s1BC', 's2BC', 's3BC')] <- sBCScores[[1]]
#   data.temp$siBC <- sBCScores[[2]]
#   data.temp[,c('s1D','s2D','s3D')] <- sDegreeScores[[1]]
#   data.temp$siD <- sDegreeScores[[2]]
#   data.temp[,c('s1HC', 's2HC', 's3HC')] <- sHCScores[[1]]
#   data.temp$siHC <- sHCScores[[2]]
#   pairwiseNetwork <- graph_from_data_frame(d = pairList, directed = FALSE, vertices = popData$ID[which(popData$Alive == "Y")])
#   E(pairwiseNetwork)$weight <- 1
#   pairwiseNetwork <- simplify(pairwiseNetwork, edge.attr.comb="sum")
#   data.temp$strength <- as.vector(strength(pairwiseNetwork))
#   data.temp$degree <- as.vector(degree(pairwiseNetwork))
#   pairwiseNetwork_InvertedWeights <- graph_from_adjacency_matrix(invert_weighted_graph(graph = pairwiseNetwork), 
#                                                                  mode = "undirected", weighted = TRUE)
#   data.temp$betweenness <- as.vector(betweenness(pairwiseNetwork_InvertedWeights, directed = FALSE, normalized = TRUE, weights = edge.attributes(pairwiseNetwork_InvertedWeights)$weight))
#   data.temp$harmonicCent <- as.vector(harmonic_centrality(pairwiseNetwork_InvertedWeights, normalized = TRUE, weights = edge.attributes(pairwiseNetwork_InvertedWeights)$weight))
# }
#   else{
#     data.temp <- data.frame(
#                           "ID" = popData$ID[which(popData$Alive == "Y")], 
#                           "timeStep" = t, 
#                           "Age" = popData$Age[which(popData$Alive == "Y")],
#                           "Parent" = popData$Parent[which(popData$Alive == "Y")],
#                           "GSPref" = popData$GSPref[which(popData$Alive == "Y")],
#                           "AgeBias" = popData$AgeBias[which(popData$Alive == "Y")],
#                           "NVariants" = get_number_of_variants_known(popData = popData, popKnowledge = popKnowledge),
#                           "NSimpleVariants" = get_number_of_variants_known(popData = popData, popKnowledge = popKnowledge, complexity = 1),
#                           "N2Variants" = get_number_of_variants_known(popData = popData, popKnowledge = popKnowledge, complexity = 2),
#                           "variantMax" = length(fullKnowledge),
#                           "simpleVariantMax" = sum(sapply(seq(1:length(fullKnowledge)), function(x) sum(sapply(fullKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == 1)))),
#                           "C2VariantMax" = sum(sapply(seq(1:length(fullKnowledge)), function(x) sum(sapply(fullKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == 2)))),
#                           "s1BC" = 0,
#                           "s2BC" = 0,
#                           "s3BC" = 0,
#                           "s4BC" = 0,
#                           "s5BC" = 0,
#                           "s6BC" = 0,
#                           "s7BC" = 0,
#                           "s8BC" = 0,
#                           "s9BC" = 0,
#                           "s10BC" = 0,
#                           "siBC" = 0,
#                           "strength" = 0, 
#                           "wBC" = 0, 
#                           "s1D" = 0,
#                           "s2D" = 0,
#                           "s3D" = 0,
#                           "s4D" = 0,
#                           "s5D" = 0,
#                           "s6D" = 0,
#                           "s7D" = 0,
#                           "s8D" = 0,
#                           "s9D" = 0,
#                           "s10D" = 0,
#                           "siD" = 0, 
#                           "s1HC" = 0,
#                           "s2HC" = 0,
#                           "s3HC" = 0,
#                           "s4HC" = 0,
#                           "s5HC" = 0,
#                           "s6HC" = 0,
#                           "s7HC" = 0,
#                           "s8HC" = 0,
#                           "s9HC" = 0,
#                           "s10HC" = 0,
#                           "siHC" = 0)
#   dualHyperNetwork <- get_dual_hypergraph(hypNet = currentPartners, popData = popData)
#   sBCScores <- get_s_betweenness(hypergraph = dualHyperNetwork, smax = 10)
#   sDegreeScores <- get_s_degree(hypergraph = dualHyperNetwork, smax = 10)
#   sHCScores <- get_s_harmonic_centrality(hypergraph = dualHyperNetwork, smax = 10)
#   weightedGraph <- get_weighted_graph(hypNet = dualHyperNetwork)
#   #binaryGraph <- convert_weighted_to_binary_graph(graph = weightedGraph)
#   invertedGraph <- graph_from_adjacency_matrix(invert_weighted_graph(graph = weightedGraph), mode = "undirected", weighted = TRUE)
#   data.temp$wBC <- as.vector(betweenness(invertedGraph, directed = FALSE, normalized = TRUE, weights = edge.attributes(invertedGraph)$weight))
#   data.temp$strength <- as.vector(strength(graph_from_adjacency_matrix(weightedGraph, mode = "undirected", weighted = TRUE), mode = "all"))
#   data.temp[,c('s1BC', 's2BC', 's3BC', 's4BC', 's5BC', 's6BC', 's7BC', 's8BC', 's9BC', 's10BC')] <- sBCScores[[1]]
#   data.temp$siBC <- sBCScores[[2]]
#   data.temp[,c('s1D','s2D','s3D','s4D','s5D','s6D','s7D','s8D','s9D','s10D')] <- sDegreeScores[[1]]
#   data.temp$siD <- sDegreeScores[[2]]
#   data.temp[,c('s1HC', 's2HC', 's3HC', 's4HC', 's5HC', 's6HC', 's7HC', 's8HC', 's9HC', 's10HC')] <- sHCScores[[1]]
#   data.temp$siHC <- sHCScores[[2]]
#   }
#   return(data.temp)
# }

convert_weighted_to_binary_graph <- function(graph) {
  matrixTemp <- graph
  matrixTemp <- sapply(seq(1:nrow(matrixTemp)), 
                       function(x) sapply(seq(1:ncol(matrixTemp)), 
                                          function(y) ifelse(matrixTemp[x,y] > 0, 1, 0)))
  rownames(matrixTemp) <- rownames(graph)
  colnames(matrixTemp) <- colnames(graph)
  return(matrixTemp)
}

invert_weighted_graph <- function(graph) {
  matrixTemp <- as_adjacency_matrix(graph, type = "both", attr = "weight", sparse = FALSE)
  maxWeight <- max(E(graph)$weight)
  matrixTemp <- sapply(seq(1:nrow(matrixTemp)), 
                       function(x) sapply(seq(1:ncol(matrixTemp)), 
                                          function(y) ifelse(matrixTemp[x,y] > 0, (1 + maxWeight) - matrixTemp[x,y], 0)))
  rownames(matrixTemp) <- rownames(graph)
  colnames(matrixTemp) <- colnames(graph)
  return(matrixTemp)
}

# record_global_summary_data <- function(popData, popKnowledge, timeStep) {
#   livePopKnowledge <- unlist(popKnowledge[popData$ID[which(popData$Alive == "Y")]])
#   uniquePopKnowledge <- unique(livePopKnowledge)
#   knowledgeRichness <- length(uniquePopKnowledge)
#   diversityVect <- sapply(seq(1:length(uniquePopKnowledge)), function(x) sum(livePopKnowledge == uniquePopKnowledge[x])/length(livePopKnowledge))
#   knowledgeDiversity <- -sum(diversityVect * log(diversityVect))
#   knowledgeEvenness <- knowledgeDiversity/log(length(uniquePopKnowledge))
#   return(c(timeStep, knowledgeRichness, knowledgeDiversity, knowledgeEvenness))
# }

# simple_contagion <- function(hyperEdges, popKnowl, popData, variantSelection = c("proportional"), probTrans) {
#   livePop <- popData[which(popData$Alive == "Y"),]
#   newKnowledge <- vector("list", nrow(livePop))
#   names(newKnowledge) <- c(livePop$ID)
#   variantVect <- sapply(seq(1:length(hyperEdges)), function(x) 
#     ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0,
#       sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1), "NA"))
#   variantReps <- sapply(seq(1:length(variantVect)), function(x) 
#     ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0, 
#            sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]),0))
#   learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
#     variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
#     0, ifelse(variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] == "NA", 0, 1)))
#   for(i in 1:length(learningStatus)) {
#     nDemons <- variantReps[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
#     newKnowledge[as.character(livePop$ID[i])] <- ifelse(learningStatus[i] == 1 & 
#                                                           sample(c(0,1), 1, 
#                                                                  prob = c((1 - probTrans) ^ nDemons, 
#                                                                           1 - (1 - probTrans)^nDemons)) == 1, 
#                                                         variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))], 
#                                                         list())
#   }
#   return(newKnowledge)
# }
# 
# complex_contagion_fractional_threshold <- function(hyperEdges, popKnowl, popData, variantSelection = c("proportional"), threshold) {
#   livePop <- popData[which(popData$Alive == "Y"),]
#   newKnowledge <- vector("list", nrow(livePop))
#   names(newKnowledge) <- c(livePop$ID)
#   variantVect <- sapply(seq(1:length(hyperEdges)), function(x) 
#     ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0,
#            sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1), "NA"))
#   variantReps <- sapply(seq(1:length(variantVect)), function(x) 
#     ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0, 
#            sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]),0))
#   learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
#     variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
#     0, ifelse(variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] == "NA", 0, 1)))
#   for(i in 1:length(learningStatus)) {
#     nDemons <- variantReps[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
#     hyperEdgeSize <- length(hyperEdges[[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]])
#     newKnowledge[as.character(livePop$ID[i])] <- ifelse(learningStatus[i] == 1 & 
#                                                           nDemons/hyperEdgeSize >= threshold, 
#                                                         variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))], 
#                                                         list())
#   }
#   return(newKnowledge)
# }
# 
# complex_contagion_numeric_threshold <- function(hyperEdges, popKnowl, popData, variantSelection = c("proportional"), threshold) {
#   livePop <- popData[which(popData$Alive == "Y"),]
#   newKnowledge <- vector("list", nrow(livePop))
#   names(newKnowledge) <- c(livePop$ID)
#   variantVect <- sapply(seq(1:length(hyperEdges)), function(x) 
#     ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0,
#            sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1), "NA"))
#   variantReps <- sapply(seq(1:length(variantVect)), function(x) 
#     ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0, 
#            sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]),0))
#   learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
#     variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
#     0, ifelse(variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] == "NA", 0, 1)))
#   for(i in 1:length(learningStatus)) {
#     nDemons <- variantReps[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
#     newKnowledge[as.character(livePop$ID[i])] <- ifelse(learningStatus[i] == 1 & 
#                                                           nDemons >= threshold, 
#                                                         variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))], 
#                                                         list())
#   }
#   return(newKnowledge)
# }

get_s_degree <- function(hypergraph, smax, vertexNames, mode) {
  matrixTemp <- matrix(0, nrow = ncol(hypergraph), ncol = smax)
  matrixTemp <- sapply(seq_along(1:smax), function(x) 
    matrixTemp[,x] <- as.vector(
      degree(
        graph_from_adjacency_matrix(
          get_s_line_graph(
            hypNet = hypergraph, size = x, vertexNames = vertexNames, mode = mode), 
          mode = "undirected"))))
  rankedMatrix <- matrix(0, nrow = ncol(hypergraph), ncol= smax)
  rankedMatrix <- sapply(seq_along(1:smax), function(x) rankedMatrix[,x] <- 
                           rank(matrixTemp[,x]))
  integratedScores <- rank(-rowSums(rankedMatrix)/smax, ties.method = "average")
  return(list(matrixTemp, integratedScores))
}

#Should add options for directed and weighted incidence matrices
get_incidence_matrix <- function(hyperNetwork, vertices) {
  hyperNetwork_Sorted <- lapply(hyperNetwork, sort)
  incidenceMatrix <- matrix(0, nrow = length(vertices), ncol = length(unique(hyperNetwork_Sorted)))
  row.names(incidenceMatrix) <- vertices
  colnames(incidenceMatrix) <- 1:length(unique(hyperNetwork_Sorted))
  for(i in 1:length(unique(hyperNetwork_Sorted))) {
    hyperEdge <- unique(hyperNetwork_Sorted)[[i]]
    hyperEdgeInstances <- sum(sapply(seq(from = 1, to = length(hyperNetwork_Sorted)), function(x) isTRUE(all.equal(hyperEdge, hyperNetwork_Sorted[[x]]))))
    for(j in hyperEdge) {
      incidenceMatrix[as.character(j), as.character(i)] <- hyperEdgeInstances
    }
  }
  return(incidenceMatrix)
}

get_dual_incidence_matrix <- function(hyperNetwork, vertices) {
  hyperNetwork_Sorted <- lapply(hyperNetwork, sort)
  incidenceMatrix <- matrix(0, nrow = length(vertices), ncol = length(hyperNetwork_Sorted))
  row.names(incidenceMatrix) <- vertices
  colnames(incidenceMatrix) <- as.integer(names(hyperNetwork_Sorted))
  for(i in 1:length(hyperNetwork_Sorted)) {
    hyperEdge <- hyperNetwork_Sorted[[i]]
    hyperEdgeInstances <- sum(sapply(seq(from = 1, to = length(hyperNetwork_Sorted)), function(x) isTRUE(all.equal(hyperEdge, hyperNetwork_Sorted[[x]]))))
    for(j in hyperEdge) {
      incidenceMatrix[as.character(j), as.character(i)] <- hyperEdgeInstances
    }
  }
  return(incidenceMatrix)
}

get_adjacency_matrix_from_incidence_matrix <- function(I, V, weighted = NULL) {
  if(weighted == TRUE) {
    binaryI <- ifelse(I[,] > 0, 1, 0)
    adjMat <- binaryI %*% 
      diag(colSums(I)/sapply(seq(from = 1, to = ncol(I)), function(x) sum(I[,x] > 0))) %*%
      t(binaryI) - diag(rowSums(I))
  } else{
    #Add check to see if hyperedges are weighted and produce warning if so
    if(sum(unique(as.vector(I)) > 1) > 0) print("Warning: hyperedges appear to be weighted. Set weighted = TRUE")
    adjMat <- I %*% t(I) - diag(rowSums(I))
  }
  row.names(adjMat) <- V
  colnames(adjMat) <- V
  return(adjMat)
}

do_social_contagion_initial <- function(incidMat, demonID, popData, 
                                lambda, v, method) {
  hyperEdgesTemp <- as.vector(which(incidMat[as.character(demonID),]>0))
  hyperEdgeWeights <- as.vector(incidMat[as.character(demonID),hyperEdgesTemp])
  if(length(hyperEdgesTemp) > 1) {
    selectedEdge <- sample(hyperEdgesTemp, size = 1, prob = hyperEdgeWeights)
  } else{
    selectedEdge <- hyperEdgesTemp
  }
  edgeMembers <- as.integer(names(which(incidMat[,selectedEdge]>0)))
  naiveEdgeMembers <- edgeMembers[edgeMembers %in% popData$ID[which(popData$knowledgeState == 0)]]
  if(length(naiveEdgeMembers) > 0) {
    if(method == "groupSizeIndependent") {
      learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < 1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v)
      learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
    } else{
      if(method == "groupSizeDependent") {
        learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < (1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v))/(length(edgeMembers) - 1)
        learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
      }
    }
  } else{
    learners <- naiveEdgeMembers
  }
  return(list(learners, selectedEdge))
}

#This is an alternative version that selects hyperedges for contagion based on their intersection with the active hyperedges on the previous time step
#For the Phil. Trans. manuscript, we only used a hyperedge weight-based selection method (do_social_contagion_initial)
# do_social_contagion <- function(incidMat, demonID, popData, informedData,
#                                         lambda, v, method) {
#   #Find membership from hyperedge in previous time step
#   edgeMembers <- as.integer(names(which(incidMat[,informedList[informedList$informedIDs == demonID, "selectedEdge"]] > 0)))
#   #Find intersecting hyperedges for demonstrator k
#   adjHyperEdges <- as.vector(which(incidMat[as.character(demonID),] > 0))
#   if(length(edgeMembers) > 1) {
#   hyperEdgeIntersections <- as.vector(colSums(incidMat[as.character(edgeMembers),] > 0))
#   adjIntersections <- hyperEdgeIntersections[adjHyperEdges]
#   originalEdgeIndex <- which(adjHyperEdges == informedList[informedList$informedIDs == demonID, "selectedEdge"])
#   adjHyperEdges <- adjHyperEdges[-originalEdgeIndex]
#   adjIntersections <- adjIntersections[-originalEdgeIndex]
#   hyperEdgeWeights <- as.vector(incidMat[as.character(demonID),adjHyperEdges])
#   
#   #This does not account for hyperedge weights, only intersections; could potentially combine by multiplying intersection size with hyperedge weight
#   if(length(adjHyperEdges) > 1) {
#       selectedEdge <- sample(adjHyperEdges, size = 1, prob = adjIntersections)
#   } else{
#     selectedEdge <- adjHyperEdges
#   }
#   } else {
#     adjIntersections <- rep(1, length(adjHyperEdges) - 1)
#     originalEdgeIndex <- which(adjHyperEdges == informedList[informedList$informedIDs == demonID, "selectedEdge"])
#     adjHyperEdges <- adjHyperEdges[-originalEdgeIndex]  #This does not account for hyperedge weights, only intersections; could potentially combine by multiplying intersection size with hyperedge weight
#     if(length(adjHyperEdges) > 1) {
#         selectedEdge <- sample(adjHyperEdges, size = 1)
#     } else{
#       selectedEdge <- adjHyperEdges
#     }
#   }
#   
#   edgeMembers <- as.integer(names(which(incidMat[,selectedEdge]>0)))
#   naiveEdgeMembers <- edgeMembers[edgeMembers %in% popData$ID[which(popData$knowledgeState == 0)]]
#   if(length(naiveEdgeMembers) > 0) {
#     if(method == "groupSizeIndependent") {
#       learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < 1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v)
#       learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
#     } else{
#       if(method == "groupSizeDependent") {
#         learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < (1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v))/(length(edgeMembers) - 1)
#         learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
#       }
#     }
#   } else{
#     learners <- naiveEdgeMembers
#   }
#   return(list(learners, selectedEdge))
# }

get_subedge_density <- function(hypergraph, focalEdge) {
  hyperEdgeMembership <- as.integer(names(which(hypergraph[,focalEdge] > 0)))
  if(length(hyperEdgeMembership) > 3) {
    possibleSubsets <- unlist(sapply(seq(from = 2, to = (length(hyperEdgeMembership) - 1)), function (x)
      combn(hyperEdgeMembership, m = x, simplify = FALSE)), recursive = FALSE)
  } else{
    possibleSubsets <- combn(hyperEdgeMembership, m = length(hyperEdgeMembership) - 1, simplify = FALSE)
  }
  hyperEdgeSizes <- sapply(seq(from = 1, to = ncol(hypergraph)), function(x) sum(hypergraph[,x] > 0))
  indicesTemp <- which(hyperEdgeSizes < length(hyperEdgeMembership) & hyperEdgeSizes > 1)
  subSetsPresent <- rep(0, length(indicesTemp))
  for(j in 1:length(indicesTemp)) {
    tempMemb <- as.integer(names(which(hypergraph[,indicesTemp[j]] > 0)))
    tempSubs <- Filter(function(x) length(x) == length(tempMemb), possibleSubsets)
    subSetsPresent[j] <- sum(sapply(seq(from = 1, to = length(tempSubs)), function(x) all(sort(tempMemb) %in% sort(tempSubs[[x]]))))
  }
  sed <- (sum(subSetsPresent)+1)/(length(possibleSubsets) + 1)
  return(sed)
}

get_global_subedge_density <- function(hypergraph) {
  hyperEdgeSizes <- sapply(seq(from = 1, to = ncol(hypergraph)), function(x) sum(hypergraph[,x] > 0))
  higherOrderIndices <- which(hyperEdgeSizes >= 3)
  subEdgeDensities <- rep(0, length(higherOrderIndices))
  for(i in 1:length(higherOrderIndices)) {
    subEdgeDensities[i] <- get_subedge_density(hypergraph = hypergraph, focalEdge = higherOrderIndices[i])
  }
  gsed <- mean(subEdgeDensities)
  return(gsed)
}

get_local_subedge_density <- function(hypergraph, vertex) {
  hyperEdgeSizes <- sapply(seq(from = 1, to = ncol(hypergraph)), function(x) sum(hypergraph[,x] > 0))
  edgeParticipation <- which(sapply(seq(from = 1, to = ncol(hypergraph)), function(x) 
    vertex %in% as.integer(names(which(hypergraph[,x] > 0)))))
  higherOrderIndices <- edgeParticipation[which(edgeParticipation %in% which(hyperEdgeSizes >= 3))]
  subEdgeDensities <- rep(0, length(higherOrderIndices))
  if(length(higherOrderIndices) > 0){
    for(i in 1:length(higherOrderIndices)) {
      subEdgeDensities[i] <- get_subedge_density(hypergraph = hypergraph, focalEdge = higherOrderIndices[i])
    }
    gsed <- mean(subEdgeDensities)
  } else{
    gsed <- 0
  }
  return(gsed)
}

select_friends <- function(prefMatrix, popData, t) {
  #livePop <- popData[which(popData$Alive == "Y"),]
  livePop <- popData
  groupList <- vector("list", nrow(livePop))
  pairList <- matrix(0, nrow = nrow(livePop), ncol = 2)
  selectOrder <- sample(livePop$ID, size = nrow(livePop), replace = FALSE)
  prefMatrixTemp <- prefMatrix[[t]]
  for(i in as.vector(which(rowSums(prefMatrixTemp) == 0))) {
    prefMatrixTemp[i,] <- 1
    prefMatrixTemp[i, i] <- 0
  }
  for(i in selectOrder) {
    if(length(groupList[sapply(groupList, function(x) i %in% x)]) == 0) {
      #currentAge <- livePop$Age[which(livePop$ID == i)]
      selectivity <- 1 + 1/(1 + exp(-(livePop$Selectivity[which(livePop$ID == i)]) * 50 + 5))
      partnerProb <- (prefMatrixTemp[as.character(i),]^selectivity)/sum((prefMatrixTemp[as.character(i),]^selectivity))
      partner <- sample(livePop$ID, size = 1, replace = FALSE, prob = partnerProb)

        if(length(groupList[sapply(groupList, function(x) partner %in% x)]) == 0) {
          groupList[[1 + (length(groupList) - length(groupList[which(sapply(groupList, is.null))]))]] <- c(i, partner)
          pairList[which(livePop$ID == i),] <- c(i,partner)
        } else{
            groupList[[which(sapply(groupList, function(x) partner %in% x))[1]]] <- c(groupList[[which(sapply(groupList, function(x) partner %in% x))[1]]], i)
            pairList[which(livePop$ID == i),] <- c(i,partner)
        }
    }
  }
  groupList <- groupList[-which(sapply(groupList, is.null))]
  pairList <- pairList[pairList[,1] > 0, ]
  return(list(groupList, pairList))
}

update_friendship_matrix <- function(prefMatrices, popData, timeStep, currentPartners) {
  currentMatrix <- prefMatrices[[timeStep]]
  for(i in 1:nrow(currentPartners)) {
    partner1 <- currentPartners[i,1]
    partner2 <- currentPartners[i,2]
    currentMatrix[as.character(partner1), as.character(partner2)] <-
      currentMatrix[as.character(partner1), as.character(partner2)] + 1
    currentMatrix[as.character(partner2), as.character(partner1)] <-
      currentMatrix[as.character(partner2), as.character(partner1)] + 1
  }
  # if(length(demoChanges[[2]]) > 0) {
  #   maxID <- max(popData$ID[which(popData$Alive == "Y")])
  #   newIDs <- rep(0, length(demoChanges[[2]]))
  #   for(i in 1:length(demoChanges[[2]])) {
  #     newIDs[i] <- maxID + i
  #   }
  #   oldNames <- colnames(currentMatrix)
  #   newIDColumns <- matrix(0, nrow = nrow(currentMatrix), ncol = length(newIDs))
  #   colnames(newIDColumns) <- as.character(newIDs)
  #   rownames(newIDColumns) <- row.names(currentMatrix)
  #   currentMatrix <- cbind(currentMatrix, newIDColumns)
  #   newIDRows <- matrix(0, nrow = length(newIDs), ncol = ncol(currentMatrix))
  #   rownames(newIDRows) <- as.character(newIDs)
  #   colnames(newIDRows) <- colnames(currentMatrix)
  #   currentMatrix <- rbind(currentMatrix, newIDRows)
  # 
  #   if(inheritance == "parental") {
  #     for(i in 1:length(newIDs)) {
  #       parent <- demoChanges[[2]][[i]]
  #       currentEdges <- unique(c(currentPartners[which(currentPartners[,1] == parent),],
  #                                currentPartners[which(currentPartners[,2] == parent),]))
  #       partners <- currentEdges[which(currentEdges != parent)]
  #       currentMatrix[as.character(newIDs[i]), ] <- 1
  #       currentMatrix[,as.character(newIDs[i])] <- 1
  #       for(j in partners) {
  #         currentMatrix[as.character(newIDs[i]), as.character(j)] <- 3
  #         currentMatrix[as.character(j), as.character(newIDs[i])] <- 3
  #       }
  #       currentMatrix[as.character(newIDs[i]), as.character(parent)] <- 5
  #       currentMatrix[as.character(parent), as.character(newIDs[i])] <- 5
  #     }
  #   }
  # 
  #   for(i in demoChanges[[1]]) {
  #     currentMatrix <- currentMatrix[!rownames(currentMatrix) %in% as.character(i), !colnames(currentMatrix) %in% as.character(i)]
  #   }
  # }
  diag(currentMatrix) <- 0
  return(currentMatrix)
}

push <- function(Q, element) {
  Q <- append(Q, list(element))
  # #Q <- Q[order(sapply(Q, function(x) x[[1]]))]
  # #Believe that I need to add the as.numeric for this to function as intended
  #Q <- Q[order(as.numeric(sapply(Q, function(x) x[[1]])))]
  Q <- Q[order(as.numeric(map(Q, 1)))]
  # Q <- rbindlist(list(Q,data.table("priorityQueue" = list(element))))
  # Q <- Q[order(as.numeric(map(Q[,priorityQueue], 1)))]
  return(Q)
}

pop <- function(Q) {
  element <- Q[[1]]
  Q <- Q[-1]
  #return(list(Q,list(as.numeric(element[[1]]), element[[2]], element[[3]], element[[4]])))
  return(list(Q,c(element)))
  # element <- Q[1,]
  # Q <- Q[-1,]
  # return(list(Q, element))
}



dijkstra <- function(G, s) {
  
  #S will contain a sequence of vertices in non-decreasing distance from source node s
  #S <- vector("character")
  
  #Suppress warning when creating data.table
  # # oldw <- getOption("warn")
  # # options(warn = -1)
  # 
  # #Create data table that includes for each vertex ID (ID)
  # #its immediate predecessors in all shortest paths from source node s (P)
  # #the distance of the shortest path from the source node s to that vertex (D)
  # #the number of shortest paths from the source node s to the vertex (sigma)
  # #and whether a vertext has been visited by the algorithm yet
  # # djkTable <- data.table(ID = sort(names(G)),
  # #                        P = list(),
  # #                        D = numeric(),
  # #                        sigma = 0,
  # #                        Visited = 0)
  
  IDList <- sort(names(G))
  lengthIDList <- length(IDList)
  P <- vector("list", lengthIDList)
  names(P) <- IDList
  djkMatrix <- matrix(c(rep(NA,lengthIDList),rep(0,2*lengthIDList)), nrow = lengthIDList, ncol = 3)
  rownames(djkMatrix) <- IDList
  colnames(djkMatrix) <- c("D", "sgm", "V")
  #djkMatrix
  # djkTable <- data.table(ID = unique(G[order(ID1),ID1]),
  #                        P = list(), 
  #                        D = numeric(),
  #                        sigma = 0,
  #                        Visited = 0)
  # on.exit(options(warn = oldw))
  
  #Initialize by setting the number of shortest paths for the source node s equal to 1
  #djkTable[ID==s,sigma := 1]
  djkMatrix[s,2] <- 1
  
  #Q will be our priority queue
  Q <- list()
  #Q <- data.table("priorityQueue" = list())

  #seen will indicate the distance of the currently known shortest path from s to each other node
  seen <- list()
  seen[[s]] <- 0
  
  #To be honest, upon reviewing the algorithm... I'm not certain what next_c is doing and whether it is needed
  #may be a holdover from the conversion from python
  c <- 1
  
  next_c <- function() {
    c <<- c + 1
    return(c - 1)
  }
  
  #I wonder if I can modify this to use data.table, as I think it will be much simpler with a data.frame
  
  #Populate the priority queue with the source node s
  #Each entry includes the distance to the source node, the focal predecessor of the node and the vertext name
  Q <- push(Q, c(as.numeric(0), as.numeric(next_c()), s, s))
  
  #Repeat so long as the priority queue is not empty
  while (length(Q) > 0) {
  # while(nrow(Q) > 0) {
    
    #Pop off the top element of the priority queue (sorted by shortest distance)
    poppedQ <- pop(Q)
    
    #Update priority queue without popped element
    Q <- poppedQ[[1]]
    #element <- flatten(poppedQ[[2]][,priorityQueue])
    
    #From the popped element, extract distance (of its focal predecessor?), the focal predecessor, and vertex name
    element <- poppedQ[[2]]
    dist <- as.numeric(element[1])
    pred <- element[3]
    v <- element[4]
    
    #If the node under consideration has already been visited by the algorithm, skip it
    # if (djkTable[ID == v, Visited] == 1) {
    if(djkMatrix[v,3] == 1) {
      next
    }
    
    #update for vertex v its current sigma value, equal to the number of shortest paths associated with v so far
    #plus all shortest paths associated with the focal predecessor of v
    
    # focalRow <- which(djkTable[,ID == v])
    # newSigma <- djkTable[ID == v, sigma] + djkTable[ID == pred, sigma]
    # set(djkTable, focalRow, "sigma", newSigma)
    #djkTable[ID == v, sigma := djkTable[ID == v, sigma] +  djkTable[ID == pred, sigma]]
    djkMatrix[v,2] <- djkMatrix[v,2] + djkMatrix[pred,2]
    
    #Add v to the list of vertices in non-decreasing distance from source node s
    #S <- c(S, v)
    
    #Update current distance of vertex v from source node s
    #djkTable[ID == v, D := as.numeric(dist)]
    #set(djkTable, focalRow, "D", as.numeric(dist))
    djkMatrix[v,1] <- dist
    #Indicate that vertex v has been visited
    #djkTable[ID == v, Visited := 1]
    #set(djkTable, focalRow, "Visited", 1)
    djkMatrix[v,3] <- 1

    # for (w in names(G[[v]])) {
    #   vw_dist <- dist + G[[v]][[w]]
    #   if (djkTable[ID == w, Visited] == 0 && (!(w %in% names(seen)) || vw_dist < seen[[w]])) {
    #     seen[[w]] <- vw_dist
    #     Q <- push(Q, c(round(vw_dist, 15), next_c(), v, w))
    #     djkTable[ID == w, sigma := 0]
    #     djkTable[ID == w, P := list(c(v))]
    #   } else if (vw_dist == seen[[w]]) {
    #     djkTable[ID == w, sigma := djkTable[ID == w, sigma] + djkTable[ID == v, sigma]]
    #     djkTable[ID == w, P := unlist(c(djkTable[ID == w, P], v))]
    #   }
    # }
    
    #Extract each node v is connected to, along with the distances from node v, and add to these the current shortest distance to v (?)  
    vw_dist_vect <- map2_dbl(pluck(G,v), dist, sum)
    # vw_dist_vect <- G[ID1 == v, Weight] + dist
    # names(vw_dist_vect) <- G[ID1 == v, ID2]
    
    #For each node connected to v
    for(w in 1:length(vw_dist_vect)){
      
      wName <- names(vw_dist_vect)[w]
      wDist <- as.vector(vw_dist_vect[w])
      # wRow <- which(djkTable[,ID == wName])
      
      #If a node has not yet been visited by the algorithm and either its name is not in 'seen' or the distance to w from source node s via v is less than the currently known shortest path to w
      #if (djkTable[ID == names(vw_dist_vect)[w], Visited] == 0 && (!(names(vw_dist_vect)[w] %in% names(seen)) || vw_dist_vect[w] < seen[[w]])) {
      #if (djkTable[ID == wName, Visited] == 0 && (!(wName %in% names(seen)) || vw_dist_vect[w] < seen[[w]])) {  
      if (djkMatrix[wName,3] == 0 && (!(wName %in% names(seen)) || wDist < seen[[w]])) {
        
        #Update the shortest path distance associated with w in seen
        # seen[[wName]] <- vw_dist_vect[w]
        seen[[wName]] <- wDist
        
        #Push into the priority queue an entry containing the shortest distance from s to w, , its focal predecessor (i.e., v), and the name of w
        Q <- push(Q, c(as.numeric(round(wDist, 15)), as.numeric(next_c()), v, wName))
        
        #Initialize the number of shortest paths containing node w to 0
        #djkTable[ID == names(vw_dist_vect)[w], sigma := 0]
        #set(djkTable,wRow,"sigma",0)
        
        #Add v to the list of predecessors of w
        #djkTable[ID == names(vw_dist_vect)[w], P := list(c(v))]
        #set(djkTable,wRow,"P",list(c(v)))
        P[[wName]] <- c(v)
        
        #If the distance captured by vw_dist_vect equals the current known shortest path to w:
      #} else if (vw_dist_vect[w] == seen[[wName]]) {
      } else if (wDist == seen[[wName]]) {
        
        #Add to the number of shortest paths associated with w all the shortest paths to its predecessor v
        #djkTable[ID == names(vw_dist_vect)[w], sigma := djkTable[ID == names(vw_dist_vect)[w], sigma] + djkTable[ID == v, sigma]]
        #sigmaTemp <- djkTable[ID == wName, sigma] + djkTable[ID == v, sigma]
        #set(djkTable,wRow,"sigma",sigmaTemp)
        djkMatrix[wName,2] <- djkMatrix[wName,2] + djkMatrix[v,2]
        
        #Add node v to the list of predecessors of w
        #djkTable[ID == names(vw_dist_vect)[w], P := unlist(c(djkTable[ID == names(vw_dist_vect)[w], P], v))]
        # predTemp <- unlist(c(djkTable[ID == wName, P],v))
        # set(djkTable,wRow,"P",predTemp)
        P[[wName]] <- c(unlist(P[wName]),v)
      }
    }
  }
  
  # The initialization of the algorithm makes sigma values double.
  # So we can return the exact values by dividing by 2.
  #djkTable[,sigma := djkTable[,sigma]/2]
  # sigmaList <- djkTable[,sigma]/2
  # set(djkTable,j = "sigma", value = sigmaList)
  djkMatrix[,2] <- djkMatrix[,2]/2
  
  # djkTable <- data.table(ID = unique(G[order(ID1),ID1]),
  #                        P = list(), 
  #                        D = numeric(),
  #                        sigma = 0,
  #                        Visited = 0)
  
  djkTable <- data.table(ID = IDList, P = P, D = djkMatrix[,1], sigma = djkMatrix[,2], Visited = djkMatrix[,3])
  return(djkTable[order(D)])
}

function4 <- function(node, window, nodes, graph, timestamp, alpha) {
  betweenness <- numeric(length(nodes))
  names(betweenness) <- sort(nodes)
  
  dummyNode <- paste(node, -1, sep = ".")
  nodes <- c(nodes, dummyNode)
  dummyEdges <- data.table(ID1 = character(), ID2 = character(), time1 = integer(), time2 = integer(), Weight = numeric())
  for(t in (timestamp - window + 1):timestamp){
    if(paste(node,t,sep = ".") %in% nodes) {
      dummyEdge <- data.table(ID1 = node, ID2 = node, time1 = -1, time2 = t, Weight = 0)
      dummyEdges <- rbindlist(list(dummyEdge, dummyEdges))
      #dummyEdges <- rbind(dummyEdge, dummyEdges)
    }
  }
  E_prime_u <- rbindlist(list(graph, dummyEdges))
  
  node_prime <- dummyNode
  
  E_prime_u$ID1.t <- paste(E_prime_u$ID1, E_prime_u$time1, sep = ".")
  E_prime_u$ID2.t <- paste(E_prime_u$ID2, E_prime_u$time2, sep = ".")

  # E_prime_u_same <- E_prime_u[time1 == time2]
  # E_prime_u_diff <- E_prime_u[time1 != time2]
  # G <- vector("list", length(unique(E_prime_u[,c(ID1.t,ID2.t)])))
  # names(G) <- unique(E_prime_u[,c(ID1.t,ID2.t)])
  # for(r in names(G)) {
  #   subNames <- c(E_prime_u_same[ID1.t == r, ID2.t], 
  #                 E_prime_u_same[ID2.t == r, ID1.t],
  #                 E_prime_u_diff[ID1.t == r, ID2.t])
  #   subWeights <- c(E_prime_u_same[ID1.t == r, Weight], 
  #                   E_prime_u_same[ID2.t == r, Weight],
  #                   E_prime_u_diff[ID1.t == r, Weight])
  #   G[[r]] <- flatten(list(subWeights))
  #   names(G[[r]]) <- subNames
  # }
  G<-list()
  for(r in 1:nrow(E_prime_u)) {
    ID1 <- E_prime_u[r,"ID1.t"]
    ID2 <- E_prime_u[r,"ID2.t"]
    if(E_prime_u$time1[r] == E_prime_u$time2[r]) {
      #G[[ID1]][ID2] <- list(E_prime_u[r,"Weight"])
      G[[unlist(ID1)]][unlist(ID2)] <- E_prime_u[r,"Weight"]
      G[[unlist(ID2)]][unlist(ID1)] <- E_prime_u[r,"Weight"]
    } else {
      G[[unlist(ID1)]][unlist(ID2)] <- E_prime_u[r,"Weight"]
      #Does this work? Update: Why is this here...??
      #G[[ID2]][ID2] <- 1
    }
  }
  # G_table <- data.table(ID1 = character(), ID2 = character(), Weight = numeric())
  # for(r in 1:nrow(E_prime_u)) {
  #   ID1 <- unlist(E_prime_u[r,"ID1.t"])
  #   ID2 <- unlist(E_prime_u[r,"ID2.t"])
  #   if(E_prime_u$time1[r] == E_prime_u$time2[r]) {
  #     G_temp <- data.table(ID1 = c(ID1,ID2), ID2 = c(ID2,ID1), Weight = unlist(E_prime_u[r,"Weight"]))
  #   } else {
  #     G_temp <- data.table(ID1 = ID1, ID2 = ID2, Weight = unlist(E_prime_u[r,"Weight"]))
  #   }
  #   G_table <- bind_rows(G_table,G_temp)
  # }
  
  result <- dijkstra(G=G, s=node_prime)
  S <- result[!(is.na(D)),ID]
  P <- result[!(is.na(D)),P]
  names(P) <- S
  #Sigma indicates number of shortest paths
  sigma <- result[!(is.na(D)),sigma]
  names(sigma) <- S
  D <- result[!(is.na(D)),D]
  
  D_prime <- numeric()
  S_prime <- character(0)
  #Sigma_prime indicates number of shortest-fastest paths
  sigma_prime <- numeric(length(S))
  names(sigma_prime) <- S
  
  for (x in S) {
    nodeTemp <- sub("\\..*", "", x)
    #Had to modify the following from *..
    timeTemp <- as.numeric(sub('.*\\.', "", x))
    if (timeTemp != -1 && (!(nodeTemp %in% names(D_prime)) || result[ID==x,D] == D_prime[nodeTemp])) {
      D_prime[nodeTemp] <- result[ID==x,D]
      S_prime <- c(S_prime, x)
      sigma_prime[x] <- result[ID==x,sigma]
    } else {
      sigma_prime[x] <- 0
    }
  }
  
  betweenness <- brandes_algo(betweenness = betweenness, S = S, P = P, sigma = sigma, sigma_prime = sigma_prime, s = node)
  return(betweenness)
}

brandes_algo <- function(betweenness,S, P, sigma, sigma_prime,s) {
  delta <- numeric(length(sigma))
  names(delta) <- names(sigma)
  
  while(length(S) > 0){
    x <- S[length(S)]
    S <- S[-length(S)]
    coeff <- ((sigma_prime[x]/sigma[x])+delta[x])/sigma[x]
    nodeTemp <- sub("\\..*", "", x)
    timeTemp <- as.numeric(sub('.*\\.', "", x))
    for(v in P[[x]]) {
      pNodeTemp <- sub("\\..*", "", v)
      pTimeTemp <- as.numeric(sub('.*\\.', "", v))
      if(pTimeTemp != -1) {
        delta[v] <- delta[v] + (sigma[v] * coeff)
      }
    }
    if(nodeTemp != s) {
      betweenness[x] <- betweenness[x] + delta[x]
    }
  }
  return(betweenness)
}

#Probably should call this function something else (get_temporal_transformation perhaps?)
get_temporal_edge_list <- function(edgeList, window, timestamp, alpha = NULL) {
  W = seq(from = timestamp - (window - 1), to = timestamp)
  
  E_prime <- data.frame("ID1" = character(),
                        "ID2" = character(),
                        "Time1" = numeric(),
                        "Time2" = numeric(),
                        "Weight" = numeric())
  
  # Update E' by taking the union of pairs ((u, t), (v, t), alpha)
  E_prime <- edgeList[which(edgeList$time1 >= min(W) & edgeList$time2 <= timestamp),]
  E_prime$Weight <- alpha
  
  vertexTuples <- sort(unique(c(paste(E_prime$ID1,E_prime$time1,sep = ","),paste(E_prime$ID2,E_prime$time1,sep = ","))))
  ID = sub(",.*", "", vertexTuples)
  Time <- as.numeric(sub(".*,", "", vertexTuples))
  vertexCopies <- data.frame("ID" = ID, "Time" = Time)
  vertexCopies <- vertexCopies[order(vertexCopies$ID, vertexCopies$Time),]
  
  # Further update E' by taking the union of pairs ((v, t), (v, t'), (1 - alpha)(t' - t))
  E_self_IDs <- as.vector(unlist(sapply(unique(vertexCopies$ID), function(x) rep(x,nrow(vertexCopies[ID == x,])-1))))
  
  E_self_time1 <- as.vector(unlist(sapply(unique(E_self_IDs), function(x) vertexCopies[ID == x,"Time"][1:(length(vertexCopies[ID == x,"Time"])-1)])))
  E_self_time2 <- as.vector(unlist(sapply(unique(E_self_IDs), function(x) vertexCopies[ID == x,"Time"][2:length(vertexCopies[ID == x,"Time"])])))
  E_self <- data.frame("ID1" = E_self_IDs, "ID2" = E_self_IDs, 
                       "time1" = E_self_time1, "time2" = E_self_time2)
  E_self$Weight <- (1-alpha)*(E_self$time2 - E_self$time1)
  
  E_prime <- rbind(E_prime, E_self)
  return(E_prime)
}

get_temporal_betweenness_centrality <- function(edgeList, alpha, window, timestamp, normalized = FALSE, method, nCores = NULL){
  
  E_prime <- get_temporal_edge_list(edgeList = edgeList, window = window, timestamp = timestamp, alpha = alpha)
  vertices <- sort(unique(c(E_prime$ID1, E_prime$ID2)))
  TBCList <- vector("list", length(vertices))
  names(TBCList) <- vertices
  
  nodes = sort(unique(c(paste(E_prime$ID1, E_prime$time1, sep = "."), paste(E_prime$ID2, E_prime$time2, sep = "."))))
  
  if(method == "parallel") {
    cl <- makeCluster(nCores)
    registerDoParallel(cores = nCores)
    #clusterEvalQ(cl)
    #clusterExport(cl, envir = .GlobalEnv, varlist = c('dijkstra', 'function4', 'brandes_algo'))
    TBCList <- foreach(i = vertices, .packages = c("data.table")) %dopar% {
      
      function4 <- function(node, window, nodes, graph, timestamp, alpha) {
        betweenness <- numeric(length(nodes))
        names(betweenness) <- sort(nodes)
        
        dummyNode <- paste(node, -1, sep = ".")
        nodes <- c(nodes, dummyNode)
        dummyEdges <- data.table("ID1" = character(), "ID2" = character(), "time1" = numeric(), "time2" = numeric(), "Weight" = numeric())
        for(t in (timestamp - window + 1):timestamp){
          if(paste(node,t,sep = ".") %in% nodes) {
            dummyEdge <- data.table("ID1" = node, "ID2" = node, "time1" = -1, "time2" = t, "Weight" = 0)
            dummyEdges <- rbind(dummyEdge, dummyEdges)
          }
        }
        E_prime_u <- rbind(graph, dummyEdges)
        
        node_prime <- dummyNode
        
        E_prime_u$ID1.t <- paste(E_prime_u$ID1, E_prime_u$time1, sep = ".")
        E_prime_u$ID2.t <- paste(E_prime_u$ID2, E_prime_u$time2, sep = ".")
        G <- list()
        for(r in 1:nrow(E_prime_u)) {
          ID1 <- E_prime_u[r,"ID1.t"]
          ID2 <- E_prime_u[r,"ID2.t"]
          if(E_prime_u$time1[r] == E_prime_u$time2[r]) {
            G[[ID1]][ID2] <- list(E_prime_u[r,"Weight"])
            G[[ID2]][ID1] <- list(E_prime_u[r,"Weight"])
          } else {
            G[[ID1]][ID2] <- list(E_prime_u[r,"Weight"])
            #Does this work? Update: Why is this here...??
            #G[[ID2]][ID2] <- 1
          }
        }
        
        result <- dijkstra(G=G, s=node_prime)
        S <- result[!(is.na(D)),ID]
        P <- result[!(is.na(D)),P]
        names(P) <- S
        #Sigma indicates number of shortest paths
        sigma <- result[!(is.na(D)),sigma]
        names(sigma) <- S
        D <- result[!(is.na(D)),D]
        
        D_prime <- numeric()
        S_prime <- character(0)
        #Sigma_prime indicates number of shortest-fastest paths
        sigma_prime <- numeric(length(S))
        names(sigma_prime) <- S
        
        for (x in S) {
          nodeTemp <- sub("\\..*", "", x)
          #Had to modify the following from *..
          timeTemp <- as.numeric(sub('.*\\.', "", x))
          if (timeTemp != -1 && (!(nodeTemp %in% names(D_prime)) || result[ID==x,D] == D_prime[nodeTemp])) {
            D_prime[nodeTemp] <- result[ID==x,D]
            S_prime <- c(S_prime, x)
            sigma_prime[x] <- result[ID==x,sigma]
          } else {
            sigma_prime[x] <- 0
          }
        }
        
        betweenness <- brandes_algo(betweenness = betweenness, S = S, P = P, sigma = sigma, sigma_prime = sigma_prime, s = node)
        return(betweenness)
      }
      
      brandes_algo <- function(betweenness,S, P, sigma, sigma_prime,s) {
        delta <- numeric(length(sigma))
        names(delta) <- names(sigma)
        
        while(length(S) > 0){
          x <- S[length(S)]
          S <- S[-length(S)]
          coeff <- ((sigma_prime[x]/sigma[x])+delta[x])/sigma[x]
          nodeTemp <- sub("\\..*", "", x)
          timeTemp <- as.numeric(sub('.*\\.', "", x))
          for(v in P[[x]]) {
            pNodeTemp <- sub("\\..*", "", v)
            pTimeTemp <- as.numeric(sub('.*\\.', "", v))
            if(pTimeTemp != -1) {
              delta[v] <- delta[v] + (sigma[v] * coeff)
            }
          }
          if(nodeTemp != s) {
            betweenness[x] <- betweenness[x] + delta[x]
          }
        }
        return(betweenness)
      }
      
      push <- function(Q, element) {
        Q <- append(Q, list(element))
        Q <- Q[order(sapply(Q, function(x) x[[1]]))]
        return(Q)
      }
      
      pop <- function(Q) {
        element <- Q[[1]]
        Q <- Q[-1]
        return(list(Q,list(as.numeric(element[[1]]), element[[2]], element[[3]], element[[4]])))
      }
      
      dijkstra <- function(G, s) {
        S <- vector("character")
        #Suppress warning when creating data.table
        oldw <- getOption("warn")
        options(warn = -1)
        djkTable <- data.table(ID = sort(names(G)),
                               P = list(), 
                               D = numeric(),
                               sigma = 0,
                               Visited = 0)
        on.exit(options(warn = oldw))
        djkTable[ID==s,sigma := 1]
        
        Q <- list()
        
        seen <- list()
        seen[[s]] <- 0
        c <- 1
        
        next_c <- function() {
          c <<- c + 1
          return(c - 1)
        }
        
        #I wonder if I can modify this to use data.table, as I think it will be much simpler with a data.frame
        Q <- push(Q, c(as.numeric(0), as.numeric(next_c()), s, s))
        while (length(Q) > 0) {
          poppedQ <- pop(Q)
          Q <- poppedQ[[1]]
          element <- poppedQ[[2]]
          dist <- as.numeric(element[[1]])
          pred <- element[[3]]
          v <- element[[4]]
          
          if (djkTable[ID == v, Visited] == 1) {
            next
          }
          
          djkTable[ID == v, sigma := djkTable[ID == v, sigma] +  djkTable[ID == pred, sigma]]
          S <- c(S, v)
          djkTable[ID == v, D := as.numeric(dist)]
          djkTable[ID == v, Visited := 1]
          
          for (w in names(G[[v]])) {
            vw_dist <- dist + G[[v]][[w]]
            if (djkTable[ID == w, Visited] == 0 && (!(w %in% names(seen)) || vw_dist < seen[[w]])) {
              seen[[w]] <- vw_dist
              Q <- push(Q, c(round(vw_dist, 15), next_c(), v, w))
              djkTable[ID == w, sigma := 0]
              djkTable[ID == w, P := list(c(v))]
            } else if (vw_dist == seen[[w]]) {
              djkTable[ID == w, sigma := djkTable[ID == w, sigma] + djkTable[ID == v, sigma]]
              djkTable[ID == w, P := unlist(c(djkTable[ID == w, P], v))]
            }
          }
        }
        
        # The initialization of the algorithm makes sigma values double.
        # So we can return the exact values by dividing by 2.
        djkTable[,sigma := djkTable[,sigma]/2]
        
        return(djkTable[order(D)])
      }
      
      list(function4(node = i, window = window, nodes = nodes, graph = E_prime, timestamp = timestamp, alpha = alpha))
    }
    stopCluster(cl)
    #potentially works; check if sorting vertex names is needed
    names(TBCList) <- vertices
    TBCList <- unlist(TBCList, recursive = FALSE)
  } else{
  for(i in vertices) {
    TBCList[i] <- list(function4(node = i, window = window, nodes = nodes, graph = E_prime, timestamp = timestamp, alpha = alpha))
  }
  }
  
  TBC <- numeric(length(vertices))
  names(TBC) <- vertices
  TBCMatrix <- do.call(rbind, TBCList)
  TBCTemp <- colSums(TBCMatrix)
  
  for(i in names(TBCTemp)) {
    nodeTemp <- sub("\\..*", "", i)
    TBC[nodeTemp] <- TBC[nodeTemp] + TBCTemp[i]
  }
  
  if(normalized) {
    TBC <- (TBC - min(TBC))/(max(TBC)-min(TBC))
  }
  
  return(TBC)
}
