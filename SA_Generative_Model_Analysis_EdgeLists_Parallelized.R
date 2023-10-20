source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
library(igraph)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_networkData="Sim-networkData_edgeListTopologies_paramSweep"
if(!file.exists(sim_networkData)) dir.create(sim_networkData)

popDataList <- importCSVs(path = "/home/mhasenja/scratch/SA_HyperNets/Run5/Sim-livingPopData_diffusionTopologies_paramSweep/")
edgeLists <- importCSVs(path = "/home/mhasenja/scratch/SA_HyperNets/Run5/Sim-edgeLists_diffusionTopologies_paramSweep/")

popData <- data.frame("simID" = rep(0, 100),
                      "ID" = 0,
                      "ID.SimID" = 0,
                      "timeStep" = 0,
                      "ageBias" = 0,
                      "selectGrad" = 0,
                      "GSPref" = 0, 
                      "Parent" = 0,
                      "Age" = 0,
                      "BIAge" = 0,
                      "WIAge" = 0,
                      "Degree" = 0,
                      "inDegree" = 0,
                      "outDegree" = 0,
                      "Strength" = 0,
                      "inStrength" = 0,
                      "outStrength" = 0,
                      "dirBetween" = 0,
                      "undirBetween" = 0)

foreach(i = 1:length(popDataList)) %dopar% {
  simIDTemp = ceiling(i/10)
  nameLength = nchar(names(popDataList)[[i]])
  timeTemp = as.integer(substring(names(popDataList)[[i]], nameLength - 6, nameLength - 4))
  popData$ID <- popDataList[[i]]$ID
  popData$simID <- simIDTemp
  popData$ID.SimID <- paste(popData$ID, popData$simID, sep = ".")
  popData$timeStep <- timeTemp
  popData$ageBias <- popDataList[[i]]$AgeBias
  popData$selectGrad <- popDataList[[i]]$selectGradient
  popData$GSPref <- popDataList[[i]]$GSPref
  popData$Parent <- popDataList[[i]]$Parent
  popData$Age <- popDataList[[i]]$Age
  
  currentIDs <- popDataList[[i]]$ID
  focalEdgeList <- edgeLists[[i]][,-1]
  pairMatrix <- matrix(data = 0, nrow = length(currentIDs), ncol = length(currentIDs))
  row.names(pairMatrix) <- currentIDs
  colnames(pairMatrix) <- currentIDs
  
  for(i in 1:nrow(focalEdgeList)) {
    pairMatrix[as.character(focalEdgeList[i,][1]), as.character(focalEdgeList[i,][2])] <- 1 + pairMatrix[as.character(focalEdgeList[i,][1]), as.character(focalEdgeList[i,][2])]
  }
  
  pairGraph <- graph_from_adjacency_matrix(pairMatrix, mode = "directed", diag = FALSE, weighted = TRUE)
  
  popData$Degree <- as.vector(degree(pairGraph, mode = "all"))
  popData$inDegree <- as.vector(degree(pairGraph, mode = "in"))
  popData$outDegree <- as.vector(degree(pairGraph, mode = "out"))
  popData$undirBetween <- as.vector(betweenness(pairGraph, directed = FALSE, normalized = TRUE, weights = 1/edge.attributes(pairGraph)$weight))
  popData$dirBetween <- as.vector(betweenness(pairGraph, directed = TRUE, normalized = TRUE, weights = 1/edge.attributes(pairGraph)$weight))
  popData$Strength <- as.vector(strength(pairGraph, loops = FALSE, mode = "all"))
  popData$inStrength <- as.vector(strength(pairGraph, loops = FALSE, mode = "in"))
  popData$outStrength <- as.vector(strength(pairGraph, loops = FALSE, mode = "out"))
  
  write.csv(popData, file = file.path(sim_networkData, sprintf("popData_%s_%.2f_%02i.csv", run_ID, simIDTemp, i)))
}

