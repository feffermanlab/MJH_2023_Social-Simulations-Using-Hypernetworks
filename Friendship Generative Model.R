
#Load required functions
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
library(igraph)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

#Set global parameters that do not change across simulations
#Parameters needed to generate initial hypergraph: are population size, number of hyperedges, and probability of hyperedge membership
N = 40
maxT = 100
nNetReps = 12

#Set number of simulations to run per set of parameters
numberSims = 100

#Create folders in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
incidence_mats = "Sim-incidMat_diffusionTopologies_paramSweep_tvbcTest"
sim_popData = "Sim-popData_diffusionTopologies_paramSweep_tvbcTest"
if(!file.exists(incidence_mats)) dir.create(incidence_mats)
if(!file.exists(sim_popData)) dir.create(sim_popData)

#Set random seed
set.seed(11102023)

#Set vector of random seeds, with a unique seed for each simulation
seedVect <- abs(round(rnorm((numberSims), mean = 500000, sd = 100000)))

foreach(s = 1:numberSims) %dopar% {
  set.seed(seedVect[s])
  popData <- data.frame("SimID" = s,
                        "ID" = seq(1:N),
                        "Selectivity" = runif(N, min = 0, max = 0.2))
  popNetwork <- list(c(1:N))
  
  #Create array to hold preference matrices for each time step
  preferenceMatrices <- initialize_preference_matrix(N = N, maxT = maxT, currentPartners = popNetwork)
  outputSteps <- c(rep(FALSE, maxT-5), rep(TRUE,5))
  timeStep <- 1
  
  repeat{
    #Create list to hold hypernetwork data for current time step
    currentPartners <- vector("list", 12 * nrow(popData))
    #Create empty vector to hold dyadic partner selection data
    pairList <- c()
    
    #Allow for groups to form based on individuals' dyadic preferences
    for(k in 1:nNetReps){
      partners <- select_friends(prefMatrix = preferenceMatrices, popData = popData, t = timeStep)
      indexStart <- (length(currentPartners) - length(currentPartners[which(sapply(currentPartners, is.null))])) + 1
      indexStop <- indexStart + length(partners[[1]]) - 1
      currentPartners[indexStart:indexStop] <- partners[[1]]
      pairList <- rbind(pairList, partners[[2]])
    }
    currentPartners <- currentPartners[-which(sapply(currentPartners, is.null))]
    
    preferenceMatrices[[timeStep + 1]] <- update_friendship_matrix(prefMatrices = preferenceMatrices, popData = popData, timeStep = timeStep, currentPartners = pairList)
    
    if(timeStep %in% (1:maxT)[outputSteps]){
      currentIncidenceMatrix <- get_incidence_matrix(hyperNetwork = currentPartners, vertices = popData$ID)
      write.csv(currentIncidenceMatrix, file = file.path(incidence_mats, sprintf("incidMat_%s_%02i_%03i.csv", run_ID, s, timeStep)))
      }
    
    if(timeStep == maxT) {
      break
    }
  timeStep <- timeStep + 1
  }
  
  currentIncidenceMatrix <- get_incidence_matrix(hyperNetwork = currentPartners, vertices = popData$ID)
  write.csv(currentIncidenceMatrix, file = file.path(incidence_mats, sprintf("incidMat_%s_%02i_%03i.csv", run_ID, s, timeStep)))
  write.csv(popData, file = file.path(sim_popData, sprintf("popData_%s_%02i.csv", run_ID, s)))
  
}