#Load packages and custom functions
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)

library(igraph)
library(dplyr)
library(foreach)
library(doParallel)
library(data.table)

registerDoParallel(cores = 20)

#stopifnot(dir.exists("Sim-incidMat_diffusionTopologies_paramSweep/"))

#Import incidence matrices
incidMats <- importCSVs(path = "~/scratch/IC_Monograph/SFHH/Incidence_Matrices/")
centralityData <- importCSVs(path = "~/scratch/IC_Monograph/SFHH/Centrality_Data/")
centralityData <- as.data.table(centralityData[[1]])

#Reformat incidence matrix
incidMat_C <- incidMats[[1]]
rownames(incidMat_C) <- incidMat_C[,1]
incidMat_C <- incidMat_C[,-1]
colnames(incidMat_C) <- 1:ncol(incidMat_C)

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-details_contagionResults_DirectMetrics_July11"
if(!file.exists(sim_details)) dir.create(sim_details)

#Set random seed to ensure repeatability
set.seed(07112024)

#Set simulation-wide parameters
initialInformed = c(1,3,6,10,20)
lambda = 0.025

#Set target seed strategy set and names
seedStrategySet <- c("Degree", "Strength", "sDeg_Avg", "sDeg_Avg_G", "random")

#Determines whether probability of learning decreases within increasing group/hyperedge size
groupInterferenceEffect <- c("groupSizeIndependent")

#Determines whether transmission resembles a simple or complex contagion
socialReinforcementValues <- c(1.1, 3)

#Set number of replicates per condition
r = 500

currentIDs <- rownames(incidMat_C)

#Complete round of sims for each incidence matrix
foreach(i = 1:r) %dopar% {
  
  simID = i
  
  focalData <- data.frame("simID" = simID, 
                          "ID" = currentIDs,
                          "knowledgeState" = 0,
                          "acquisitionTime" = 0, 
                          "initDemons" = 0)
  
  startingData <- focalData
  
  for (c in 1:length(socialReinforcementValues)) {
    
    for(g in 1:length(groupInterferenceEffect)) {
      
      for(s in 1:length(seedStrategySet)) {
        
        for(d in 1:length(initialInformed)) {
        
        focalData <- startingData
        seedStrategy <- seedStrategySet[[s]]
        seedSetSize <- initialInformed[[d]]
        
        if(seedStrategy != "random") {
          possibleSeeds <- centralityData[get(seedStrategy) >= as.vector(quantile(centralityData[,get(seedStrategy)],0.9))]$ID
          seedSet <- sample(possibleSeeds, seedSetSize, replace = FALSE)
        } else {
          seedSet <- sample(centralityData$ID, seedSetSize, replace = FALSE)
        }
        
        focalData[which(focalData$ID %in% seedSet),]$knowledgeState <- 1
        focalData[which(focalData$ID %in% seedSet),]$initDemons <- 1
        focalData$socialReinforcement <- socialReinforcementValues[c]
        focalData$seedStrategy <- seedStrategy
        
        v = socialReinforcementValues[c]
        t = 1
        
        informedList <- data.frame("informedIDs" = focalData$ID[which(focalData$knowledgeState == 1)])
        
        repeat{
          newLearners <- data.frame("informedIDs" = 0)
          
          for(k in informedList$informedIDs) {
            contagionOutcome <- do_social_contagion_initial(incidMat = incidMat_C, demonID = k, popData = focalData, 
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
          if(sum(focalData$knowledgeState) >= (0.95 * nrow(focalData)) | 
             t >= 5000) {
            break
          }
        }
        
        write.csv(focalData, file = file.path(sim_details, sprintf("simData_%s_%.1f_%.01i_%.01i_%.02i.csv", run_ID, v, d, s, simID)))
        
      }
      }
    }
  }
}
