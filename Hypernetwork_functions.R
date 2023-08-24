
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
    hyperEdgeIndices <- hyperEdgeIndices[which(hyperEdgeIndices != ncol(hypNet))]
    #can i modfiy this so that i and j are both the relevant hyperedgeindices
    if(length(hyperEdgeIndices) > 1) {
    for(i in hyperEdgeIndices) {
    #for(i in 1:(length(vertexNames)-1)) {
      refVect <- as.vector(which(hypNet[,i]>0))
      focalCol <- hypNet[,i]
      indices <- hyperEdgeIndices[which(hyperEdgeIndices != i)]
      matrixTemp[i,c(indices)] <- sapply(indices, 
      #matrixTemp[i,(i+1):ncol(matrixTemp)] <- sapply(seq(from = i + 1, to = ncol(matrixTemp)), 
                         function(x)
                           ifelse(
                           length(focalCol[c(intersect(refVect, 
                                                      as.vector(which(hypNet[,x]>0))))]) >= size, 1, 0))
      # for(j in (i+1):length(vertexNames)) {
      #   #if(i != j) {
      #     #Modified to allow for duplicates of the same hyperedge to count
      #     if(
      #       # length(intersect(as.vector(which(hypNet[,i]>0)), 
      #       #                   as.vector(which(hypNet[,j]>0)))) > 0 & 
      #        sum(hypNet[,i][c(intersect(as.vector(which(hypNet[,i]>0)), 
      #                                   as.vector(which(hypNet[,j]>0))))]) >= size) {
      #       matrixTemp[i,j] <- 1
      #     }
      #   #}
      # }
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
      #* associateData$RelativePrevInt
      associateData$PsAdjust <- associateData$PsUnadjust/sum(associateData$PsUnadjust)
      associateData$Ps <- associateData$PsAdjust * associateData$RelativePrevInt
      # associateData$Ps <- (popData$AgeBias[which(popData$ID == i)] * (associateData$Age/sum(associateData$Age))) + 
      #   ((1 - popData$AgeBias[which(popData$ID == i)]) * associateData$RelativePrevInt)
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

# knowledgeCapacity <- function(hypEdgeList, popData) {
#   knowledgeVect <- rep(0, length(hypEdgeList))
#   for(i in 1:length(hypEdgeList)) {
#     knowledgeVect[i] <- mean(popData$Knowl[which(popData$ID %in% hypEdgeList[[i]])])
#   }
#   return(knowledgeVect)
# }
# 
# transmit_information <- function(hypEdgeList, edgeKnowledge, popData) {
#   livePop <- popData[which(popData$Alive == "Y"),]
#   for(i in livePop$ID){
#     currentKnowledge <- livePop$Knowl[which(livePop$ID == i)]
#     hyperEdgeKnowledge <- edgeKnowledge[c(
#       sapply(hypEdgeList, function(x) i %in% x)
#     )]
#     popData$Knowl[which(popData$ID == i)] <- max(c(currentKnowledge, hyperEdgeKnowledge))
#   }
#   return(popData)
# }

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

generate_full_knowledge_set <- function(numDomains, numLevels, complexity) {
  domains <- letters[1:numDomains]
  
  #May want to modify this eventually to allow for different number of levels in each domain
  levels <- seq(from = 1, to = numLevels)
  
  simpleKnowledgeSet <- list()
  
  for(i in domains){
    for(j in levels){
      simpleKnowledgeSet <- c(simpleKnowledgeSet, 
                              paste(i, j, sep = ""))
    }
  }
  knowledgeSet <- simpleKnowledgeSet
  
  if(complexity > 1) {
    complexKnowledgeSet <- list()
    for(c in 2:complexity){
      for(i in domains){
        kTemp <- grep(i, simpleKnowledgeSet, value = TRUE)
        #complexKnowledgeSet <- c(complexKnowledgeSet, combn(kTemp, c, simplify = FALSE))
        complexKnowledgeSet <- c(complexKnowledgeSet, c(combn(kTemp, c, simplify = FALSE, FUN = paste, collapse = '')))
      }
    }
    knowledgeSet <- c(simpleKnowledgeSet, complexKnowledgeSet)
  }
  
  return(knowledgeSet)
}

generate_initial_pop_knowledge <- function(popData, knowledgeSet, asoc, simpleKnowledge = TRUE){
  popKnowledge <- vector("list", nrow(popData))
  
  #Works for simple knowledge; will need to decide how to handle complex knowledge (e.g., "a1 a2")
  if(simpleKnowledge == TRUE) {
    simpleKnowl <- knowledgeSet[sapply(knowledgeSet, function(x) length(x) == 1)]
    popAges <- popData$Age[which(popData$Alive == "Y")]
    numDiscoveries <- rep(0, length(popAges))
    for(i in 1:length(popAges)) {
      numDiscoveries[i] <- sum(sample(c(0,1), popAges[i], replace = TRUE, prob = c(1 - asoc, asoc)))
    }
    for(i in 1:length(numDiscoveries)) {
      popKnowledge[i] <- list(sample(simpleKnowl, numDiscoveries[i], replace = FALSE))
    }
  } else{
    fullKnowl <- knowledgeSet
    popAges <- popData$Age[which(popData$Alive == "Y")]
    numDiscoveries <- rep(0, length(popAges))
    for(i in 1:length(popAges)) {
      numDiscoveries[i] <- sum(sample(c(0,1), popAges[i], replace = TRUE, prob = c(1 - asoc, asoc)))
    }
    for(i in 1:length(numDiscoveries)) {
      popKnowledge[i] <- list(sample(fullKnowl, numDiscoveries[i], replace = FALSE))
    }
  }
  return(popKnowledge)
}

initialize_preference_matrix <- function(N, maxT, currentPartners) {
  #matrixList <- vector("list", maxT)
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

get_hyperedge_knowledge <- function(hyperEdges, popKnowl) {
  hypEKnowledge <- data.frame("edgeID" = 1:length(hyperEdges), 
                              "KTotal" = 0,
                              "KDivers" = 0)
  for(i in 1:length(hyperEdges)) {
    kTemp <- unlist(sapply(hyperEdges[[i]], function(x) unlist(popKnowl[[x]])))
    hypEKnowledge$KTotal[hypEKnowledge$edgeID == i] <- length(kTemp)
    hypEKnowledge$KDivers[hypEKnowledge$edgeID == i] <- length(unique(kTemp))
  }
  return(hypEKnowledge)
}

learn_new_knowledge <- function(hyperEdgeKnowledge, popData, popKnowl, hyperEdges, selectMethod, kThreshold) {
  livePop <- popData[which(popData$Alive == "Y"),]
  newKnowledge <- vector("list", nrow(livePop))
  names(newKnowledge) <- c(livePop$ID)
  if(selectMethod == "total"){
    maxTotal <- rep(0, length(livePop$ID))
  for(i in 1:length(livePop$ID)) {
    currentEdges <- c(which(sapply(hyperEdges, function (x) livePop$ID[i] %in% x)))
      maxTotal[i] <- max(hyperEdgeKnowledge$KTotal[which(hyperEdgeKnowledge$edgeID %in% currentEdges)])
  }
    for(i in which(maxTotal > 0)){
      #if(maxTotal > 0) {
        currentEdges <- c(which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x)))
        selectedEdge <- hyperEdgeKnowledge$edgeID[which(hyperEdgeKnowledge$edgeID %in% currentEdges & 
                                                          hyperEdgeKnowledge$KTotal == maxTotal[i])]
        #if(length(selectedEdge) > 1) {
          selectedEdge <- as.integer(sample(as.character(selectedEdge), size = 1))
        #}
        hyperEdgeKnowlVariants <- unlist(sapply(hyperEdges[[selectedEdge]], function(x) unlist(popKnowl[[x]])))
        variantData <- data.frame("variantID" = sort(unique(hyperEdgeKnowlVariants)), 
                                  "propKnowledgeable" = 0)
        for(r in 1:nrow(variantData)) {
          variantData[r,"propKnowledgeable"] <- sum(hyperEdgeKnowlVariants == variantData[r, "variantID"])/length(hyperEdges[[selectedEdge]])
        }
        learnableVariants <- variantData[which(variantData$propKnowledgeable > kThreshold),]
        if(nrow(learnableVariants) > 0) {
          learnedVariants <- learnableVariants$variantID[which(!(learnableVariants$variantID %in% unlist(popKnowl[[livePop$ID[i]]])))]
          if(length(learnedVariants) > 0){
            newKnowledge[as.character(livePop$ID[i])] <- list(learnedVariants)
          }
          }
      }
    } 
    else{
      if(selectMethod == "diversity") {
        uniqueTotal <- max(hyperEdgeKnowledge$KDivers[which(hyperEdgeKnowledge$edgeID %in% currentEdges)])
        selectedEdge <- hyperEdgeKnowledge$edgeID[which(hyperEdgeKnowledge$edgeID %in% currentEdges & 
                                                          hyperEdgeKnowledge$KDivers == uniqueTotal)]
        if(length(selectedEdge) > 1) {
          selectedEdge <- sample(selectedEdge, n = 1)
        } 
      }
    }
  return(newKnowledge)
}

update_knowledge <- function(popKnowl, newKnowl) {
  newKnowl[sapply(newKnowl, is.null)] <- NULL
  for(r in 1:length(newKnowl)) {
    popKnowl[[as.integer(names(newKnowl)[[r]])]] <- append(popKnowl[[as.integer(names(newKnowl)[[r]])]], as.list(newKnowl[[r]]))
  }
  return(popKnowl)
}

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
    
    #Later I could set up a model state as Graham did in the supply chain model and specify whether to use the random or parental variant without an if check
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
      #I don't believe that the following code actually links child to parent; need to rectify this
      #I believe I have modified to allow for parental connection.
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

generate_knowledge_for_new_indivs <- function(popData, knowledgeSet, asoc, simpleKnowledge = TRUE){
  knowledgeTemp <- vector("list", length(popData$ID[which(popData$Age == 1)]))
  
  #Works for simple knowledge; will need to decide how to handle complex knowledge (e.g., "a1 a2")
  if(simpleKnowledge == TRUE) {
    simpleKnowl <- knowledgeSet[sapply(knowledgeSet, function(x) length(x) == 1)]
    numDiscoveries <- sample(c(0,1), length(which(popData$Age == 1)), replace = TRUE, prob = c(1 - asoc, asoc))
    for(i in 1:length(numDiscoveries)) {
      knowledgeTemp[i] <- list(sample(simpleKnowl, numDiscoveries[i], replace = FALSE))
    }
  }
  return(knowledgeTemp)
}

asocial_learning <- function(popKnowledge, popData, knowledgeSet, simpleKnowledge, asoc) {
  if(simpleKnowledge == TRUE) {
    numDiscoveries <- sample(c(0,1), length(which(popData$Age > 1 & popData$Alive == "Y")), replace = TRUE, prob = c(1-asoc, asoc))
    livePopIDs <- popData$ID[popData$Alive == "Y" & popData$Age > 1]
    for(i in which(numDiscoveries == 1)) {
      knowledgePool <- which(!(knowledgeSet[sapply(knowledgeSet, function(x) length(x) == 1)] %in% unlist(popKnowledge[[livePopIDs[i]]])))
      if(length(knowledgePool) > 0) {
        newKnowl <- sample(knowledgePool, size = 1)
        popKnowledge[[livePopIDs[i]]] <- append(popKnowledge[[livePopIDs[i]]], knowledgeSet[[newKnowl]])
      }
    }
  }
  return(popKnowledge)
}

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

get_number_of_variants_known  <- function(popData, popKnowledge, complexity = NULL) {
  liveIDs <- popData$ID[which(popData$Alive == "Y")]
  variantData <- rep(0, length(liveIDs))
  ifelse(is.null(complexity), 
         variantData <- sapply(liveIDs, 
                               function(x) length(unlist(popKnowledge[[x]]))), 
         variantData <- sapply(liveIDs, 
                               function(x) ifelse(length(popKnowledge[[x]]) == 0, 0, 
                                                  sum(sapply(popKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == complexity))
                               )
         )
  )
  return(variantData)
}

record_individual_summary_data <- function(popData, popKnowledge = NULL, knowledgeSet = NULL, smax, mode) {
  if(is.null(knowledgeSet)) {    
    data.temp <- data.frame(
    "ID" = popData$ID[which(popData$Alive == "Y")], 
    "timeStep" = t, 
    "Age" = popData$Age[which(popData$Alive == "Y")],
    "Parent" = popData$Parent[which(popData$Alive == "Y")],
    "GSPref" = popData$GSPref[which(popData$Alive == "Y")],
    "s1BC" = 0, "s2BC" = 0, "s3BC" = 0,
    "siBC" = 0,
    "s1D" = 0, "s2D" = 0, "s3D" = 0,
    "siD" = 0, 
    "s1HC" = 0, "s2HC" = 0, "s3HC" = 0,
    "siHC" = 0, 
    "degree" = 0,
    "strength" = 0, 
    "betweenness" = 0, 
    "harmonicCent" = 0)
  dualHyperNetwork <- get_dual_hypergraph(hypNet = currentPartners, popData = popData)
  sBCScores <- get_s_betweenness(hypergraph = dualHyperNetwork, smax = smax, mode = mode)
  sDegreeScores <- get_s_degree(hypergraph = dualHyperNetwork, smax = smax, mode = mode)
  sHCScores <- get_s_harmonic_centrality(hypergraph = dualHyperNetwork, smax = smax, mode = mode)
  data.temp[,c('s1BC', 's2BC', 's3BC')] <- sBCScores[[1]]
  data.temp$siBC <- sBCScores[[2]]
  data.temp[,c('s1D','s2D','s3D')] <- sDegreeScores[[1]]
  data.temp$siD <- sDegreeScores[[2]]
  data.temp[,c('s1HC', 's2HC', 's3HC')] <- sHCScores[[1]]
  data.temp$siHC <- sHCScores[[2]]
  pairwiseNetwork <- graph_from_data_frame(d = pairList, directed = FALSE, vertices = popData$ID[which(popData$Alive == "Y")])
  E(pairwiseNetwork)$weight <- 1
  pairwiseNetwork <- simplify(pairwiseNetwork, edge.attr.comb="sum")
  data.temp$strength <- as.vector(strength(pairwiseNetwork))
  data.temp$degree <- as.vector(degree(pairwiseNetwork))
  pairwiseNetwork_InvertedWeights <- graph_from_adjacency_matrix(invert_weighted_graph(graph = pairwiseNetwork), 
                                                                 mode = "undirected", weighted = TRUE)
  data.temp$betweenness <- as.vector(betweenness(pairwiseNetwork_InvertedWeights, directed = FALSE, normalized = TRUE, weights = edge.attributes(pairwiseNetwork_InvertedWeights)$weight))
  data.temp$harmonicCent <- as.vector(harmonic_centrality(pairwiseNetwork_InvertedWeights, normalized = TRUE, weights = edge.attributes(pairwiseNetwork_InvertedWeights)$weight))
}
  else{
    data.temp <- data.frame(
                          "ID" = popData$ID[which(popData$Alive == "Y")], 
                          "timeStep" = t, 
                          "Age" = popData$Age[which(popData$Alive == "Y")],
                          "Parent" = popData$Parent[which(popData$Alive == "Y")],
                          "GSPref" = popData$GSPref[which(popData$Alive == "Y")],
                          "AgeBias" = popData$AgeBias[which(popData$Alive == "Y")],
                          "NVariants" = get_number_of_variants_known(popData = popData, popKnowledge = popKnowledge),
                          "NSimpleVariants" = get_number_of_variants_known(popData = popData, popKnowledge = popKnowledge, complexity = 1),
                          "N2Variants" = get_number_of_variants_known(popData = popData, popKnowledge = popKnowledge, complexity = 2),
                          "variantMax" = length(fullKnowledge),
                          "simpleVariantMax" = sum(sapply(seq(1:length(fullKnowledge)), function(x) sum(sapply(fullKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == 1)))),
                          "C2VariantMax" = sum(sapply(seq(1:length(fullKnowledge)), function(x) sum(sapply(fullKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == 2)))),
                          "s1BC" = 0,
                          "s2BC" = 0,
                          "s3BC" = 0,
                          "s4BC" = 0,
                          "s5BC" = 0,
                          "s6BC" = 0,
                          "s7BC" = 0,
                          "s8BC" = 0,
                          "s9BC" = 0,
                          "s10BC" = 0,
                          "siBC" = 0,
                          "strength" = 0, 
                          "wBC" = 0, 
                          "s1D" = 0,
                          "s2D" = 0,
                          "s3D" = 0,
                          "s4D" = 0,
                          "s5D" = 0,
                          "s6D" = 0,
                          "s7D" = 0,
                          "s8D" = 0,
                          "s9D" = 0,
                          "s10D" = 0,
                          "siD" = 0, 
                          "s1HC" = 0,
                          "s2HC" = 0,
                          "s3HC" = 0,
                          "s4HC" = 0,
                          "s5HC" = 0,
                          "s6HC" = 0,
                          "s7HC" = 0,
                          "s8HC" = 0,
                          "s9HC" = 0,
                          "s10HC" = 0,
                          "siHC" = 0)
  dualHyperNetwork <- get_dual_hypergraph(hypNet = currentPartners, popData = popData)
  sBCScores <- get_s_betweenness(hypergraph = dualHyperNetwork, smax = 10)
  sDegreeScores <- get_s_degree(hypergraph = dualHyperNetwork, smax = 10)
  sHCScores <- get_s_harmonic_centrality(hypergraph = dualHyperNetwork, smax = 10)
  weightedGraph <- get_weighted_graph(hypNet = dualHyperNetwork)
  #binaryGraph <- convert_weighted_to_binary_graph(graph = weightedGraph)
  invertedGraph <- graph_from_adjacency_matrix(invert_weighted_graph(graph = weightedGraph), mode = "undirected", weighted = TRUE)
  data.temp$wBC <- as.vector(betweenness(invertedGraph, directed = FALSE, normalized = TRUE, weights = edge.attributes(invertedGraph)$weight))
  data.temp$strength <- as.vector(strength(graph_from_adjacency_matrix(weightedGraph, mode = "undirected", weighted = TRUE), mode = "all"))
  data.temp[,c('s1BC', 's2BC', 's3BC', 's4BC', 's5BC', 's6BC', 's7BC', 's8BC', 's9BC', 's10BC')] <- sBCScores[[1]]
  data.temp$siBC <- sBCScores[[2]]
  data.temp[,c('s1D','s2D','s3D','s4D','s5D','s6D','s7D','s8D','s9D','s10D')] <- sDegreeScores[[1]]
  data.temp$siD <- sDegreeScores[[2]]
  data.temp[,c('s1HC', 's2HC', 's3HC', 's4HC', 's5HC', 's6HC', 's7HC', 's8HC', 's9HC', 's10HC')] <- sHCScores[[1]]
  data.temp$siHC <- sHCScores[[2]]
  }
  return(data.temp)
}

convert_weighted_to_binary_graph <- function(graph) {
  matrixTemp <- graph
  matrixTemp <- sapply(seq(1:nrow(matrixTemp)), 
                       function(x) sapply(seq(1:ncol(matrixTemp)), 
                                          function(y) ifelse(matrixTemp[x,y] > 0, 1, 0)))
  rownames(matrixTemp) <- rownames(graph)
  colnames(matrixTemp) <- colnames(graph)
  # for(i in 1:nrow(matrixTemp)){
  #   for(j in 1:ncol(matrixTemp)) {
  #     if(matrixTemp[i,j] > 0) {
  #       matrixTemp[i,j] <- 1
  #     }
  #   }
  # }
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

record_global_summary_data <- function(popData, popKnowledge, timeStep) {
  livePopKnowledge <- unlist(popKnowledge[popData$ID[which(popData$Alive == "Y")]])
  uniquePopKnowledge <- unique(livePopKnowledge)
  knowledgeRichness <- length(uniquePopKnowledge)
  diversityVect <- sapply(seq(1:length(uniquePopKnowledge)), function(x) sum(livePopKnowledge == uniquePopKnowledge[x])/length(livePopKnowledge))
  knowledgeDiversity <- -sum(diversityVect * log(diversityVect))
  knowledgeEvenness <- knowledgeDiversity/log(length(uniquePopKnowledge))
  return(c(timeStep, knowledgeRichness, knowledgeDiversity, knowledgeEvenness))
}

simple_contagion <- function(hyperEdges, popKnowl, popData, variantSelection = c("proportional"), probTrans) {
  livePop <- popData[which(popData$Alive == "Y"),]
  newKnowledge <- vector("list", nrow(livePop))
  names(newKnowledge) <- c(livePop$ID)
  variantVect <- sapply(seq(1:length(hyperEdges)), function(x) 
    ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0,
      sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1), "NA"))
  variantReps <- sapply(seq(1:length(variantVect)), function(x) 
    ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0, 
           sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]),0))
  learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
    variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
    0, ifelse(variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] == "NA", 0, 1)))
  for(i in 1:length(learningStatus)) {
    #variantTemp <- variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
    nDemons <- variantReps[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
    newKnowledge[as.character(livePop$ID[i])] <- ifelse(learningStatus[i] == 1 & 
                                                          sample(c(0,1), 1, 
                                                                 prob = c((1 - probTrans) ^ nDemons, 
                                                                          1 - (1 - probTrans)^nDemons)) == 1, 
                                                        variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))], 
                                                        list())
  }
  return(newKnowledge)
}

complex_contagion_fractional_threshold <- function(hyperEdges, popKnowl, popData, variantSelection = c("proportional"), threshold) {
  livePop <- popData[which(popData$Alive == "Y"),]
  newKnowledge <- vector("list", nrow(livePop))
  names(newKnowledge) <- c(livePop$ID)
  variantVect <- sapply(seq(1:length(hyperEdges)), function(x) 
    ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0,
           sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1), "NA"))
  variantReps <- sapply(seq(1:length(variantVect)), function(x) 
    ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0, 
           sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]),0))
  learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
    variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
    0, ifelse(variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] == "NA", 0, 1)))
  # variantVect <- sapply(seq(1:length(hyperEdges)), function(x) sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1))
  # variantReps <- sapply(seq(1:length(variantVect)), function(x) sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]))
  # learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
  #   variantVect[which(sapply(hyperEdges, function(y) x %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
  #   0,1))
  for(i in 1:length(learningStatus)) {
    nDemons <- variantReps[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
    hyperEdgeSize <- length(hyperEdges[[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]])
    newKnowledge[as.character(livePop$ID[i])] <- ifelse(learningStatus[i] == 1 & 
                                                          nDemons/hyperEdgeSize >= threshold, 
                                                        variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))], 
                                                        list())
  }
  return(newKnowledge)
}

complex_contagion_numeric_threshold <- function(hyperEdges, popKnowl, popData, variantSelection = c("proportional"), threshold) {
  livePop <- popData[which(popData$Alive == "Y"),]
  newKnowledge <- vector("list", nrow(livePop))
  names(newKnowledge) <- c(livePop$ID)
  variantVect <- sapply(seq(1:length(hyperEdges)), function(x) 
    ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0,
           sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1), "NA"))
  variantReps <- sapply(seq(1:length(variantVect)), function(x) 
    ifelse(sum(as.vector(sapply(hyperEdges[[x]], function(y) length(popKnowl[[y]])))) > 0, 
           sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]),0))
  learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
    variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
    0, ifelse(variantVect[which(sapply(hyperEdges, function(y) livePop$ID[x] %in% y))] == "NA", 0, 1)))
  # variantVect <- sapply(seq(1:length(hyperEdges)), function(x) sample(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))), 1))
  # variantReps <- sapply(seq(1:length(variantVect)), function(x) sum(as.vector(sapply(hyperEdges[x], function(y) unlist(popKnowl[y]))) == variantVect[x]))
  # learningStatus <- sapply(seq(1:length(livePop$ID)), function(x) ifelse(
  #   variantVect[which(sapply(hyperEdges, function(y) x %in% y))] %in% unlist(popKnowl[livePop$ID[x]]),
  #   0,1))
  for(i in 1:length(learningStatus)) {
    nDemons <- variantReps[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))]
    newKnowledge[as.character(livePop$ID[i])] <- ifelse(learningStatus[i] == 1 & 
                                                          nDemons >= threshold, 
                                                        variantVect[which(sapply(hyperEdges, function(x) livePop$ID[i] %in% x))], 
                                                        list())
  }
  return(newKnowledge)
}



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

get_s_degreeV2 <- function(hypergraph, smax, vertexNames, mode) {
  matrixTemp <- matrix(0, nrow = 
                         100
                       #length(hypergraph)
                       , ncol = smax)
  matrixTemp <- sapply(seq_along(1:smax), function(x) 
    matrixTemp[,x] <- as.vector(
      degree(
        graph_from_adjacency_matrix(
          get_s_line_graph(
            hypergraph, size = x, vertexNames = vertexNames, mode = mode), 
          mode = "undirected"))))
  rankedMatrix <- matrix(0, nrow = 100
                         #length(hypergraph)
                         , ncol= smax)
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
                                #newLearners, 
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
      # if(length(learners) > 0){
      #   newLearners <- append(newLearners, learners)
      # }
    } else{
      if(method == "groupSizeDependent") {
        learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < (1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v))/(length(edgeMembers) - 1)
        learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
        # if(length(learners) > 0){
        #   newLearners <- append(newLearners, learners)
        # } 
      }
    }
  } else{
    learners <- naiveEdgeMembers
  }
  return(list(learners, selectedEdge))
}


do_social_contagion <- function(incidMat, demonID, popData, informedData,
                                        #newLearners, 
                                        lambda, v, method) {
  #Find membership from hyperedge in previous time step
  edgeMembers <- as.integer(names(which(incidMat[,informedList[informedList$informedIDs == demonID, "selectedEdge"]] > 0)))
  #Find intersecting hyperedges for demonstrator k
  adjHyperEdges <- as.vector(which(incidMat[as.character(demonID),] > 0))
    #as.vector(which(colSums(incidMat[as.character(edgeMembers),] > 0) > 0))
  if(length(edgeMembers) > 1) {
  hyperEdgeIntersections <- as.vector(colSums(incidMat[as.character(edgeMembers),] > 0))
  adjIntersections <- hyperEdgeIntersections[adjHyperEdges]
  originalEdgeIndex <- which(adjHyperEdges == informedList[informedList$informedIDs == demonID, "selectedEdge"])
  adjHyperEdges <- adjHyperEdges[-originalEdgeIndex]
  adjIntersections <- adjIntersections[-originalEdgeIndex]
  #hyperEdgesTemp <- as.vector(which(incidMat[as.character(demonID),]>0))
  hyperEdgeWeights <- as.vector(incidMat[as.character(demonID),adjHyperEdges])
  
  #This does not account for hyperedge weights, only intersections; could potentially combine by multiplying intersection size with hyperedge weight
  if(length(adjHyperEdges) > 1) {
      selectedEdge <- sample(adjHyperEdges, size = 1, prob = adjIntersections)
  } else{
    selectedEdge <- adjHyperEdges
  }
  } else {
    #hyperEdgeIntersections <- rep(1, length(adjHyperEdges))
    adjIntersections <- rep(1, length(adjHyperEdges) - 1)
    originalEdgeIndex <- which(adjHyperEdges == informedList[informedList$informedIDs == demonID, "selectedEdge"])
    adjHyperEdges <- adjHyperEdges[-originalEdgeIndex]  #This does not account for hyperedge weights, only intersections; could potentially combine by multiplying intersection size with hyperedge weight
    if(length(adjHyperEdges) > 1) {
        selectedEdge <- sample(adjHyperEdges, size = 1)
    } else{
      selectedEdge <- adjHyperEdges
    }
  }
  
  # adjIntersections <- hyperEdgeIntersections[adjHyperEdges]
  # originalEdgeIndex <- which(adjHyperEdges == informedList[informedList$informedIDs == demonID, "selectedEdge"])
  # adjHyperEdges <- adjHyperEdges[-originalEdgeIndex]
  # adjIntersections <- adjIntersections[-originalEdgeIndex]
  # #hyperEdgesTemp <- as.vector(which(incidMat[as.character(demonID),]>0))
  # hyperEdgeWeights <- as.vector(incidMat[as.character(demonID),adjHyperEdges])
  # 
  # #This does not account for hyperedge weights, only intersections; could potentially combine by multiplying intersection size with hyperedge weight
  # if(length(adjHyperEdges > 1)) {
  #   if(length(adjIntersections) == sum(adjIntersections))
  #   selectedEdge <- sample(adjHyperEdges, size = 1, prob = adjIntersections)
  # } else{
  #   selectedEdge <- adjHyperEdges
  # }
  edgeMembers <- as.integer(names(which(incidMat[,selectedEdge]>0)))
  naiveEdgeMembers <- edgeMembers[edgeMembers %in% popData$ID[which(popData$knowledgeState == 0)]]
  if(length(naiveEdgeMembers) > 0) {
    if(method == "groupSizeIndependent") {
      learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < 1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v)
      learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
      # if(length(learners) > 0){
      #   newLearners <- append(newLearners, learners)
      # }
    } else{
      if(method == "groupSizeDependent") {
        learnAttempts <- runif(n = length(naiveEdgeMembers), min = 0, max = 1) < (1 - (1 - lambda)^((length(edgeMembers) - length(naiveEdgeMembers))^v))/(length(edgeMembers) - 1)
        learners <- naiveEdgeMembers[which(learnAttempts == TRUE)]
        # if(length(learners) > 0){
        #   newLearners <- append(newLearners, learners)
        # } 
      }
    }
  } else{
    learners <- naiveEdgeMembers
  }
  return(list(learners, selectedEdge))
}

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