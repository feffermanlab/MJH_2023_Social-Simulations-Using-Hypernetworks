#Load packages and custom functions
source("Hypernetwork_functions.R")
source("dataImportAndAggregation.R")
source("seedStrategyFunctions.R")
library(igraph)
library(dplyr)

#Set working directory to where incidence matrices are stored
setwd("C://Users/matth/Desktop/Betweenness in Hypernetworks/Social ageing model - Diffusion topologies/Sim-incidMat_diffusionTopologies/")

#Import incidence matrices
incidMats <- importCSVs()

setwd("C://Users/matth/Desktop/Betweenness in Hypernetworks/Social ageing model - Diffusion topologies/Sim-edgeLists_diffusionTopologies/")
edgeLists <- importCSVs()

#Set working directory to where simulation results should be stored
setwd("C://Users/matth/Desktop/Higher-order contagion model - Simulation results/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-details_higherOrderContagion_groupSizeIndependent_Weights_PairvsHyp"
if(!file.exists(sim_details)) dir.create(sim_details)

#Set random seed to ensure repeatbility
set.seed(4122023)

#Set simulation-wide parameters
initialInformed = 3
lambda = 0.025

seedStrategySet <- c(seedStrategy_highestDegree,seedStrategy_highestBetweenness,
                     seedStrategy_highestsiD, seedStrategy_highestsiBC,
                     seedStrategy_highestGroupsiD,seedStrategy_highestGroupsiBC)

seedStrategyNames <- c("highestDegree", "highestBetweenness", 
                       "highestsiD", "highestsiBC", 
                       "highestGroupsiD", "highestGroupsiBC")

socialReinforcementValues <- c(1, 3)

#Complete round of sims for each incidence matrix
for(i in 1:length(incidMats)) {
  
  #Set up appropriate (hyper)graphs for current incidence matrix
  currentIDs <- incidMats[[i]][,1]
  focalIncidMat <- as.matrix(incidMats[[i]][,-1])
  row.names(focalIncidMat) <- currentIDs
  colnames(focalIncidMat) <- seq(from = 1, to = dim(focalIncidMat)[2])
  dualIncidMat <- t(focalIncidMat)
  GoGmatrix <- get_adjacency_matrix_from_incidence_matrix(I = focalIncidMat, 
                                                          V = currentIDs, 
                                                          weighted = TRUE)
  GoGgraph <- graph_from_adjacency_matrix(GoGmatrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  pairGraph <- graph.data.frame(d = edgeLists[[i]][,-1], directed = FALSE, vertices = currentIDs)
  E(pairGraph)$weight <- 1
  pairGraph <- simplify(pairGraph, edge.attr.comb = list(weight="sum"))
  
  simID = i
  
  focalData <- data.frame("simID" = simID, 
                          "ID" = currentIDs,
                          "degree" = 0,
                          "betweenness" = 0, 
                          "siD" = 0,
                          "siBC" = 0,
                          "topGroupsiD" = 0,
                          "topGroupsiBC" = 0, 
                          "knowledgeState" = 0,
                          "acquisitionTime" = 0, 
                          "initDemons" = 0)
  focalData$degree <- as.vector(degree(pairGraph))
  focalData$betweenness <- as.vector(betweenness(pairGraph, 
                                                 directed = FALSE, normalized = TRUE, 
                                                 weights = 1/edge.attributes(pairGraph)$weight))
  focalData$siD <- get_s_degree(hypergraph = dualIncidMat, smax = 4, vertexNames = currentIDs, mode = "incidence")[[2]]
  focalData$siBC <- get_s_betweenness(hypergraph = dualIncidMat, smax = 4, vertexNames = currentIDs, mode = "incidence")[[2]]
  groupRanks_siD <- get_s_degree(hypergraph = focalIncidMat, smax = 4, vertexNames = seq(1:dim(focalIncidMat)[2]), mode = "incidence")[[2]]
  focalData$topGroupsiD <- which(groupRanks_siD == min(groupRanks_siD))[1]
  groupRanks_siBC <- get_s_betweenness(hypergraph = focalIncidMat, smax = 4, vertexNames = seq(1:dim(focalIncidMat)[2]), mode = "incidence")[[2]]
  focalData$topGroupsiBC <- which(groupRanks_siBC == min(groupRanks_siBC))[1]
  
  startingData <- focalData
  
  for (c in 1:length(socialReinforcementValues)) {
    
    for(s in 1:length(seedStrategySet)) {
      
      focalData <- startingData
      seedStrategyFunction = seedStrategySet[[s]]
      focalData <- seedStrategyFunction(data = focalData, numSeeds = initialInformed)
      focalData$socialReinforcement <- socialReinforcementValues[c]
      focalData$seedStrategy <- seedStrategyNames[s]
      
      v = socialReinforcementValues[c]
      t = 1
      
      informedList <- data.frame("informedIDs" = focalData$ID[which(focalData$knowledgeState == 1)])
      
      repeat{
        newLearners <- data.frame("informedIDs" = 0)
        
        for(k in informedList$informedIDs) {
          contagionOutcome <- do_social_contagion_initial(incidMat = focalIncidMat, demonID = k, popData = focalData, 
                                                          lambda = lambda, v = v, method = "groupSizeIndependent")
          
          if(length(unlist(contagionOutcome[[1]])) > 0) {
            learnerTemp <- data.frame("informedIDs" = unlist(contagionOutcome[[1]]))
            newLearners <- rbind(newLearners, learnerTemp)
          }
        }
        
        newLearners <- newLearners[-1,,drop = FALSE]
        
        if(nrow(newLearners) > 0) {
          newLearners <- newLearners[!duplicated(newLearners$informedIDs),,drop = FALSE]
          focalData[which(focalData$ID %in% newLearners$informedIDs),]$knowledgeState <- 1
          focalData[which(focalData$ID %in% newLearners$informedIDs),]$acquisitionTime <- t
          informedList <- rbind(informedList, newLearners)
        }
        
        t <- t + 1
        if(sum(focalData$knowledgeState) == max(components(GoGgraph)$csize) | 
           t >= 5000){
          break
        }
      }
      
      write.csv(focalData, file = file.path(sim_details, sprintf("simData_%s_%.1f_%.01i_%.02i.csv", run_ID, v, s, simID)))
      
    }
  }
}


