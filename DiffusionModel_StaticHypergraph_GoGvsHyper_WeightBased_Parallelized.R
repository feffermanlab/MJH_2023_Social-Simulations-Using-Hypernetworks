#Load packages and custom functions
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/seedStrategyFunctions.R", chdir = TRUE)
library(igraph)
library(dplyr)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

#stopifnot(dir.exists("Sim-incidMat_diffusionTopologies_paramSweep/"))

#Import incidence matrices
incidMats <- importCSVs(path = "~/scratch/SA_HyperNets/Run7/Sim-incidMat_diffusionTopologies_paramSweep/")
popData <- importCSVs(path = "~/scratch/SA_HyperNets/Run7/Sim-livingPopData_diffusionTopologies_paramSweep/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-details_higherOrderContagion_Weights_GoGvsHyp_NetMets_Mar24"
if(!file.exists(sim_details)) dir.create(sim_details)

#Set random seed to ensure repeatability
set.seed(10092023)

#Set simulation-wide parameters
initialInformed = 3
lambda = 0.025

#Set target seed strategy set and names

#seedStrategySet <- c(seedStrategy_highestDegree, seedStrategy_highestBetweenness, seedStrategy_highestStrength,
#  seedStrategy_highestsiD, seedStrategy_highestsiBC, seedStrategy_highestSED)
#seedStrategySet <- c(seedStrategy_highestDegree, seedStrategy_highestStrength, seedStrategy_highestsiD)
#seedStrategySet <- c(seedStrategy_Oldest, seedStrategy_randomAge)
seedStrategySet <- c(seedStrategy_highestStrength, seedStrategy_highestsiD)

#seedStrategyNames <- c("highestDegree", "highestBetweenness", "highestStrength", "highestsiD", "highestsiBC", "highestSED")
#seedStrategyNames <- c("highestDegree", "highestStrength", "highestsiD")
#seedStrategyNames <- c("Oldest", "randomAge")
seedStrategyNames <- c("highestStrength", "highestsiD")

#Determines whether probability of learning decreases within increasing group/hyperedge size
groupInterferenceEffect <- c("groupSizeIndependent", "groupSizeDependent")

#Determines whether transmission resembles a simple or complex contagion
socialReinforcementValues <- c(1.1, 3)

#Complete round of sims for each incidence matrix
foreach(i = 1:length(incidMats)) %dopar% {
  
  #Set up appropriate (hyper)graphs for current incidence matrix
  currentIDs <- incidMats[[i]][,1]
  #IDMatch <- sum(currentIDs == popData[[i]]$ID)
  focalIncidMat <- as.matrix(incidMats[[i]][,-1])
  #focalIncidMat <- as.matrix(focalIncidMat[,-dim(focalIncidMat)[2]])
  row.names(focalIncidMat) <- currentIDs
  colnames(focalIncidMat) <- seq(from = 1, to = dim(focalIncidMat)[2])
  dualIncidMat <- t(focalIncidMat)
  GoGmatrix <- get_adjacency_matrix_from_incidence_matrix(I = focalIncidMat, 
                                                          V = currentIDs, 
                                                          weighted = TRUE)
  GoGgraph <- graph_from_adjacency_matrix(GoGmatrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  simID = i
  
  focalData <- data.frame("simID" = simID, 
                          "ageBias" = as.numeric(substring(names(incidMats[i]), 27, 30)),
                          "selectGrad" = as.numeric(substring(names(incidMats[i]), 32, 35)),
                          "ID" = currentIDs,
                          "degree" = 0,
                          "betweenness" = 0, 
                          "strength" = 0,
                          "siD" = 0,
                          "siBC" = 0,
                          "subEdgeDens" = 0,
                          "knowledgeState" = 0,
                          "acquisitionTime" = 0, 
                          "initDemons" = 0,
                          #"IDMatch" = IDMatch,
                          "GSPref" = popData[[i]]$GSPref)

  focalData$degree <- as.vector(degree(GoGgraph))
  focalData$betweenness <- as.vector(betweenness(GoGgraph, directed = FALSE, normalized = TRUE, weights = 1/edge.attributes(GoGgraph)$weight))
  focalData$strength <- as.vector(strength(GoGgraph, loops = FALSE, mode = "all"))
  focalData$siD <- get_s_degree(hypergraph = dualIncidMat, smax = 6, vertexNames = currentIDs, mode = "incidence")[[2]]
  focalData$siBC <- get_s_betweenness(hypergraph = dualIncidMat, smax = 6, vertexNames = currentIDs, mode = "incidence")[[2]]
  focalData$subEdgeDens <- sapply(focalData$ID, function(x) get_local_subedge_density(hypergraph = focalIncidMat, vertex = x))

  startingData <- focalData
  
  for (c in 1:length(socialReinforcementValues)) {
    
    for(g in 1:length(groupInterferenceEffect)) {
    
    for(s in 1:length(seedStrategySet)) {
      
      focalData <- startingData
      seedStrategyFunction = seedStrategySet[[s]]
      focalData <- seedStrategyFunction(data = focalData, numSeeds = initialInformed)
      focalData$socialReinforcement <- socialReinforcementValues[c]
      focalData$seedStrategy <- seedStrategyNames[s]
      focalData$groupEffect <- groupInterferenceEffect[g]
      
      v = socialReinforcementValues[c]
      t = 1
      
      informedList <- data.frame("informedIDs" = focalData$ID[which(focalData$knowledgeState == 1)])
      
      repeat{
        newLearners <- data.frame("informedIDs" = 0)
        
        for(k in informedList$informedIDs) {
          contagionOutcome <- do_social_contagion_initial(incidMat = focalIncidMat, demonID = k, popData = focalData, 
                                                          lambda = lambda, v = v, method = groupInterferenceEffect[g])
          
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
      
      write.csv(focalData, file = file.path(sim_details, sprintf("simData_%s_%.1f_%.01i_%.01i_%.02i.csv", run_ID, v, g, s, simID)))
      
    }
    }
  }
}
