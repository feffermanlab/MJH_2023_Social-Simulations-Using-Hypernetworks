#Load packages and custom functions
source("Hypernetwork_functions.R")
source("dataImportAndAggregation.R")
source("seedStrategyFunctions.R")
library(igraph)
library(dplyr)

#Set working directory to where incidence matrices are stored
setwd("C://Users/matth/Desktop/Betweenness in Hypernetworks/Social ageing model - Diffusion topologies/Sim-incidMat_diffusionTopologies/")
#setwd("C://Users/matth/Desktop/Betweenness in Hypernetworks/Social ageing model - Diffusion topologies/Sim-edgeLists_diffusionTopologies/")

#Import incidence matrices
incidMats <- importCSVs()

#Set working directory to where simulation results should be stored
setwd("C://Users/matth/Desktop/Higher-order contagion model - Simulation results/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-details_higherOrderContagion_groupSizeIndependent_Walking"
if(!file.exists(sim_details)) dir.create(sim_details)

#Set random seed to ensure repeatbility
set.seed(4122023)

#Set simulation-wide parameters
initialInformed = 3
lambda = 0.025

seedStrategySet <- c(seedStrategy_highestDegree,seedStrategy_highestBetweenness,
                     #seedStrategy_highestharmCent,
                     seedStrategy_highestsiD, seedStrategy_highestsiBC, 
                     #seedStrategy_highestsiHC,
                     seedStrategy_highestGroupsiD,seedStrategy_highestGroupsiBC
                     #seedStrategy_highestGroupsiHC
                     )

seedStrategyNames <- c("highestDegree", "highestBetweenness", 
                       #"highestHarmCent", 
                       "highestsiD", "highestsiBC", 
                       #"highestsiHC", 
                       "highestGroupsiD", "highestGroupsiBC"
                       #"highestGroupsiHC"
                       )
#contagionVariants <- c("simple", "complex")
socialReinforcementValues <- c(1, 3)

# #for sims 57 onwards
# set.seed(5423453)
# 
# #for sims 113 onwards
# set.seed(468115)

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
  # GoGgraph_InvertedWeights <- graph_from_adjacency_matrix(invert_weighted_graph(graph = GoGgraph), 
  #                                                         mode = "undirected", weighted = TRUE)
  
  simID = i
  
  focalData <- data.frame("simID" = simID, 
                          "ID" = currentIDs,
                          "degree" = 0,
                          "betweenness" = 0, 
                          #"harmCentrality" = 0, 
                          "siD" = 0,
                          "siBC" = 0,
                          #"siHC" = 0,
                          "topGroupsiD" = 0,
                          "topGroupsiBC" = 0, 
                          #"topGroupsiHC" = 0,
                          "knowledgeState" = 0,
                          "acquisitionTime" = 0, 
                          "initDemons" = 0)
  focalData$degree <- as.vector(degree(GoGgraph))
  focalData$betweenness <- as.vector(betweenness(GoGgraph, 
                                                 directed = FALSE, normalized = TRUE, 
                                                 weights = 1/edge.attributes(GoGgraph)$weight))
  # focalData$harmCentrality <- as.vector(harmonic_centrality(GoGgraph_InvertedWeights, normalized = TRUE, 
  #                                                           weights = edge.attributes(GoGgraph_InvertedWeights)$weight))
  focalData$siD <- get_s_degree(hypergraph = dualIncidMat, smax = 4, vertexNames = currentIDs, mode = "incidence")[[2]]
  focalData$siBC <- get_s_betweenness(hypergraph = dualIncidMat, smax = 4, vertexNames = currentIDs, mode = "incidence")[[2]]
  # focalData$siHC <- get_s_harmonic_centrality(hypergraph = dualIncidMat, smax = 3, vertexNames = currentIDs, mode = "incidence")[[2]]
  groupRanks_siD <- get_s_degree(hypergraph = focalIncidMat, smax = 4, vertexNames = seq(1:dim(focalIncidMat)[2]), mode = "incidence")[[2]]
  focalData$topGroupsiD <- which(groupRanks_siD == min(groupRanks_siD))[1]
  groupRanks_siBC <- get_s_betweenness(hypergraph = focalIncidMat, smax = 4, vertexNames = seq(1:dim(focalIncidMat)[2]), mode = "incidence")[[2]]
  focalData$topGroupsiBC <- which(groupRanks_siBC == min(groupRanks_siBC))[1]
  # groupRanks_siHC <- get_s_harmonic_centrality(hypergraph = focalIncidMat, smax = 3, vertexNames = seq(1:dim(focalIncidMat)[2]), mode = "incidence")[[2]]
  # focalData$topGroupsiHC <- which(groupRanks_siHC == min(groupRanks_siHC))[1]
  
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
      
      informedList <- data.frame("informedIDs" = focalData$ID[which(focalData$knowledgeState == 1)], 
                                 "selectedEdge" = 0)

        newLearners <- data.frame("informedIDs" = 0,
                                  "selectedEdge" = 0)

        #informedIDs <- focalData$ID[which(focalData$knowledgeState==1)]
        
        for(k in informedList$informedIDs) {
          contagionOutcome <- do_social_contagion_initial(incidMat = focalIncidMat, demonID = k, popData = focalData, 
                                                  #newLearners = newLearners, 
                                             lambda = lambda, v = v, method = "groupSizeIndependent")
          informedList[which(informedList$informedIDs == k),]$selectedEdge <- contagionOutcome[[2]]
          if(length(unlist(contagionOutcome[[1]])) > 0) {
            learnerTemp <- data.frame("informedIDs" = unlist(contagionOutcome[[1]]), 
                                      "selectedEdge" = contagionOutcome[[2]])
            newLearners <- rbind(newLearners, learnerTemp)
          }
        }
        
        newLearners <- newLearners[-1,]
        
        if(nrow(newLearners) > 0) {
          newLearners <- newLearners[!duplicated(newLearners$informedIDs),]
          focalData[which(focalData$ID %in% newLearners$informedIDs),]$knowledgeState <- 1
          focalData[which(focalData$ID %in% newLearners$informedIDs),]$acquisitionTime <- t
          informedList <- rbind(informedList, newLearners)
        }
        
        t <- t + 1
        
        repeat{
        #for(t in 2:11){
          newLearners <- data.frame("informedIDs" = 0,
                                    "selectedEdge" = 0)
          informedOrder <- sample(informedList$informedIDs, size = nrow(informedList), replace = FALSE)
          for(k in informedOrder) {
            # for(k in 1:71) {
            #   k = informedOrder[k]

            
            contagionOutcome <- do_social_contagion(incidMat = focalIncidMat, demonID = k, popData = focalData, informedData = informedList,
                                                            #newLearners = newLearners, 
                                                            lambda = lambda, v = v, method = "groupSizeIndependent")
            informedList[which(informedList$informedIDs == k),]$selectedEdge <- contagionOutcome[[2]]
            if(length(unlist(contagionOutcome[[1]])) > 0) {
              learnerTemp <- data.frame("informedIDs" = unlist(contagionOutcome[[1]]), 
                                        "selectedEdge" = contagionOutcome[[2]])
              newLearners <- rbind(newLearners, learnerTemp)
            }
          }
          
          newLearners <- newLearners[-1,]
          
          if(nrow(newLearners) > 0) {
            newLearners <- newLearners[!duplicated(newLearners$informedIDs),]
            focalData[which(focalData$ID %in% newLearners$informedIDs),]$knowledgeState <- 1
            focalData[which(focalData$ID %in% newLearners$informedIDs),]$acquisitionTime <- t
            informedList <- rbind(informedList, newLearners)
          }
          
        t <- t + 1
        #if(sum(focalData$knowledgeState) == nrow(focalData)) {
        if(sum(focalData$knowledgeState) == max(components(GoGgraph)$csize) | 
           t >= 5000){
          break
        }
      }
      
      write.csv(focalData, file = file.path(sim_details, sprintf("simData_%s_%.1f_%.01i_%.02i.csv", run_ID, v, s, simID)))
      
    }
  }
}


