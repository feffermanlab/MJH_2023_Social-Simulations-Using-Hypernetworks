
#Load required functions
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
library(igraph)
library(stringr)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

#Set global parameters that do not change across simulations
#Parameters needed to generate initial hypergraph: are population size, number of hyperedges, and probability of hyperedge membership
N = 100
E = 35
p = 0.05
maxT = 100
nNetReps = 12

#Set number of simulations to run per set of parameters
numberSims = 1000

#Set inheritance pattern
#Can be "random" or "parental"
#Parental means offspring are automatically connected to all individual that its parent had associated with that year
#Random means offspring are connected to a number of individuals equal to those its parent had associated with that year, randomly selected
inheritance = "parental"

#Set levels for age biases
#This governs the extent to which individuals prioritize associating with older individuals or peers
ageBiasSet = c(0, 0.5, 1.0)

#Set levels for social selectivity
#These parameters govern how rapidly an individual prioritizes established relationships over developing new ones as they age
selectiveSet = c(0, 0.10, 0.20)

#Create folders in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
edge_lists = "Sim-edgeLists_diffusionTopologies_paramSweep"
incidence_mats = "Sim-incidMat_diffusionTopologies_paramSweep"
sim_popData = "Sim-popData_diffusionTopologies_paramSweep"
livingPopData = "Sim-livingPopData_diffusionTopologies_paramSweep"
if(!file.exists(edge_lists)) dir.create(edge_lists)
if(!file.exists(incidence_mats)) dir.create(incidence_mats)
if(!file.exists(sim_popData)) dir.create(sim_popData)
if(!file.exists(livingPopData)) dir.create(livingPopData)

#ageDist <- generate_age_structure(1000, maxIter = 20000)
ageDist <- read.csv("", header = TRUE)[,2]

#Set random seed
set.seed(08262023)

#Set vector of random seeds, with a unique seed for each simulation
seedVect <- round(rnorm((numberSims * length(ageBiasSet)
                         * length(selectiveSet)
                         ), mean = 500000, sd = 200000))
seedVectList <- vector("list", length(ageBiasSet))
seedSubVectLength <- length(seedVect)/length(ageBiasSet)
for(i in 1:length(ageBiasSet)) {
  indexStart <- 1 + (i - 1)*seedSubVectLength
  indexEnd <- indexStart + seedSubVectLength - 1
  seedVectList[[i]] <- seedVect[indexStart:indexEnd]
}

foreach(s = 1:numberSims) %dopar% {

for(a in 1:length(ageBiasSet)) {
  ageBias = ageBiasSet[a]
  seedVectTemp <- as.vector(seedVectList[[a]])
  for(m in 1:length(selectiveSet)) {
    selectGradient = selectiveSet[m]

      #Set unique seed for simulation run
      set.seed(seedVectTemp[numberSims * (m - 1) + s])
      
      #Create dataframe to record data for population members
      popData <- create_population_data(N, initialAges = sample(ageDist, size = N, replace = TRUE), ageBias = ageBias, selectGradient = selectGradient)
      
      popNetwork <- generate_ER_hypergraph(v = N, m = E, p = p)
      
      #Create array to hold preference matrices for each time step
      preferenceMatrices <- initialize_preference_matrix(N = N, maxT = maxT, currentPartners = popNetwork)
      
      timeStep <- 1
      repeat{
        #Create list to hold hypernetwork data for current time step
        currentPartners <- vector("list", 12 * nrow(popData[which(popData$Alive == "Y"),]))
        #Create empty vector to hold dyadic partner selection data
        pairList <- c()
        
        #Allow for groups to form based on individuals' dyadic preferences
        for(k in 1:nNetReps){
          partners <- select_partners(prefMatrix = preferenceMatrices, popData = popData, t = timeStep)
          indexStart <- (length(currentPartners) - length(currentPartners[which(sapply(currentPartners, is.null))])) + 1
          indexStop <- indexStart + length(partners[[1]]) - 1
          currentPartners[indexStart:indexStop] <- partners[[1]]
          pairList <- rbind(pairList, partners[[2]])
        }
        currentPartners <- currentPartners[-which(sapply(currentPartners, is.null))]
        
        demoChanges <- identify_demographic_changes(popData = popData)
        
        preferenceMatrices[[timeStep + 1]] <- update_preference_matrix(prefMatrices = preferenceMatrices, popData = popData, timeStep = timeStep, 
                                                                       demoChanges = demoChanges, currentPartners = currentPartners, inheritance = inheritance)
        
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
        
        if(sum(popData$ID[which(popData$Alive == "Y")] %in% 1:N) ==0){
          break
        }
        timeStep <- timeStep + 1
      }
      
      currentPrefMat <- preferenceMatrices[[timeStep + 1]]
      preferenceMatrices <- lapply(1:maxT, matrix, data = 0, nrow = N, ncol = N)
      preferenceMatrices[[1]] <- currentPrefMat
      outputSteps <- c(rep(FALSE, maxT/2), ((maxT/2 + 1):maxT %% 5) == 0)
      
      #Begin running through time steps
      for(t in 1:maxT){
        
        #Set current time step
        timeStep <- t
        
        #Create list to hold hypernetwork data for current time step
        currentPartners <- vector("list", 12 * nrow(popData[which(popData$Alive == "Y"),]))
        #Create empty vector to hold dyadic partner selection data
        pairList <- c()
        
        #Allow for groups to form based on individuals' dyadic preferences
        for(k in 1:nNetReps){
          partners <- select_partners(prefMatrix = preferenceMatrices, popData = popData, t = timeStep)
          indexStart <- (length(currentPartners) - length(currentPartners[which(sapply(currentPartners, is.null))])) + 1
          indexStop <- indexStart + length(partners[[1]]) - 1
          currentPartners[indexStart:indexStop] <- partners[[1]]
          pairList <- rbind(pairList, partners[[2]])
          }
        currentPartners <- currentPartners[-which(sapply(currentPartners, is.null))]
        
        
        #Output pairwise edge list and hypergraph incidence matrix (only do this every 5 steps during the final 100 timesteps of simulation)
        #if(t == maxT) {
        if(t %in% (1:maxT)[outputSteps]){
          write.csv(pairList, file = file.path(edge_lists, sprintf("edgeData_%s_%.2f_%.2f_%02i_%03i.csv", run_ID, ageBias, selectGradient, s, timeStep)))
          currentIncidenceMatrix <- get_incidence_matrix(hyperNetwork = currentPartners, vertices = popData$ID[which(popData$Alive == "Y")])
          write.csv(currentIncidenceMatrix, file = file.path(incidence_mats, sprintf("incidMat_%s_%.2f_%.2f_%02i_%03i.csv", run_ID, ageBias, selectGradient, s, timeStep)))
          livingPop <- popData[which(popData$Alive == "Y"),]
          write.csv(livingPop, file = file.path(livingPopData, sprintf("livingPop_%s_%.2f_%.2f_%02i_%03i.csv", run_ID, ageBias, selectGradient, s, timeStep)))
        }
        
        #Check if any deaths occur; an individual is born for each death to keep population size stable
        demoChanges <- identify_demographic_changes(popData = popData)
        
        preferenceMatrices[[timeStep + 1]] <- update_preference_matrix(prefMatrices = preferenceMatrices, popData = popData, timeStep = timeStep, 
                                                                       demoChanges = demoChanges, currentPartners = currentPartners, inheritance = inheritance)
        
        #Update data for population members
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

      }
      
      popData <- popData[which(popData$Age > 0),]
      write.csv(popData, file = file.path(sim_popData, sprintf("popData_%s_%.2f_%.2f_%02i.csv", run_ID, ageBias, selectGradient, s)))
    }
  }
}
