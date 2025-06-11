
#Load packages
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/hypernetwork_group_size_vs_composition_simFunctions.R", chdir = TRUE)
#source("hypernetwork_group_size_vs_composition_simFunctions.R")
library(igraph)
library(asnipe)
library(data.table)
library(foreach)
library(doParallel)
library(stringr)

registerDoParallel(cores = 20)

#setwd("C://Users/matth/Desktop/Multilayer Hypernetworks/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_indData="Sim_individual_level_data"
sim_summaryData="Sim_summary_data"
sim_diffSummaryData = "Sim_diffusion_summary_data"
sim_netData_dyNet1 = "Sim_dyadic_network_1"
sim_netData_hypNet1 = "Sim_higherOrder_network_1"
sim_netData_dyNet2 = "Sim_dyadic_network_2"
sim_netData_hypNet2 = "Sim_higherOrder_network_2"
sim_netData_dyNet3 = "Sim_dyadic_network_3"
sim_netData_hypNet3 = "Sim_higherOrder_network_3"

if(!file.exists(sim_indData)) dir.create(sim_indData)
if(!file.exists(sim_summaryData)) dir.create(sim_summaryData)
if(!file.exists(sim_diffSummaryData)) dir.create(sim_diffSummaryData)
if(!file.exists(sim_netData_dyNet1)) dir.create(sim_netData_dyNet1)
if(!file.exists(sim_netData_hypNet1)) dir.create(sim_netData_hypNet1)
if(!file.exists(sim_netData_dyNet2)) dir.create(sim_netData_dyNet2)
if(!file.exists(sim_netData_hypNet2)) dir.create(sim_netData_hypNet2)
if(!file.exists(sim_netData_dyNet3)) dir.create(sim_netData_dyNet3)
if(!file.exists(sim_netData_hypNet3)) dir.create(sim_netData_hypNet3)

##Set seed for reproducibility
set.seed(5162025)

#Set social transmission rate; determines baseline probability of learning from an active neighbor in the dyadic layer
social_trans <- 0.1

#Possible radii for determinig subgroups
rad <- c(0.15, 0.25, 0.35)

#Distributions for dominance ranks
domDist <- c("domUni", "domExp")

resourceAvailability <- c(0.25, 0.5, 0.75, 1)

groupAdjust <- c("FALSE", "TRUE")

###################################################

#Set number of sims per condition
nSims = 250

foreach(s = 1:nSims) %dopar% {
  #for(s in 1:nSims) {
  
  ##Generate individuals in family/kinship groups
  #For the moment, not using families in model, but keeping in case of expanding
  ind_data <- generate_population(popSize = 100, nInitInformed = 1)
  n_indivs <- nrow(ind_data)
  ind_data$simID <- s
  
  ind_data$domUni <- runif(nrow(ind_data), min = 0, max = 1)
  domExp <- sort(rexp(nrow(ind_data), rate = 2))
  domExp <- sapply(rank(ind_data$domUni), function(x) domExp[x])
  ind_data$domExp <- (domExp - min(domExp))/(max(domExp) - min(domExp))
  ind_data$domExp <- ifelse(ind_data$domExp == 0, min(ind_data[which(ind_data$domExp > 0),]$domExp)/2, ind_data$domExp)
  
  ind_dataOrig <- ind_data
  
  ##Create spatial proximity network and hypergraphs
  netList <- generate_latent_space_multilayer_hypergraph(ind_data = ind_data, r = rad)
  
  write.csv(netList[[1]], file.path(sim_netData_dyNet1, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[2]], file.path(sim_netData_hypNet1, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[3]], file.path(sim_netData_dyNet2, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[4]], file.path(sim_netData_hypNet2, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[5]], file.path(sim_netData_dyNet3, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[6]], file.path(sim_netData_hypNet3, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  
  for(d in 1:length(domDist)) {
    for(R in resourceAvailability) {
      for(r in c(1,3,5)) {
        for(g in groupAdjust) {
        
        #Reset individual-level data to initial values
        ind_data <- ind_dataOrig
        
        domDistribution <- domDist[d]
        ind_data$distUsed <- domDistribution
        ind_data$resAvail <- R
        ind_data$network <- paste("dyadic", (r+1)/2, sep = "_")
        focaldomDist <- as.vector(unlist(ind_data[domDistribution]))
        
        resources <- round(R * mean(colSums(netList[[r+1]])))
        
        #Create list for holding output on each time step
        dataList <- vector("list", 8000)
        
        #Set initial time
        t = 1
        
        repeat{
          
          #Identify current set of informed nodes
          informedNodes <- ind_data[which(ind_data$informed == 1),]$id
          
          #Each informed individual updates its likelihood of producing a novel behavior by:
          learningOutcomes <- dyadic_diffusion(ind_data = ind_data, network = r, netList = netList, informedNodes = informedNodes, domValues = focaldomDist, groupAdjustment = g, resources = resources)
          ind_data$active <- learningOutcomes[[1]]
          ind_data$newLearners <- learningOutcomes[[2]]
          
          for(i in 1:nrow(ind_data)) {
            if(ind_data$active[i] == 1 & ind_data$firstProd[i] == 0) {
              ind_data$firstProd[i] <- t
            }
          }
          
          #Record those individuals that acquired the trait and note the time step at which this occurred
          ind_data$acqTime <- ifelse(ind_data$newLearners == 1 & ind_data$informed == 0, t, ind_data$acqTime)
          
          #Add newly informed individuals to the list of informed individuals
          ind_data$informed <- ifelse(ind_data$informed == 0, 
                                      ifelse(ind_data$newLearners == 1, 1, 0), 1)
          
          dataTemp <- data.table("simID" = s,
                                 "domDist" = domDistribution,
                                 "network" = paste("dyadic", (r+1)/2, sep = "_"),
                                 "resAvail" = R,
                                 "groupAdjust" = g,
                                 "nInformed" = 0,
                                 "percInformed" = 0,
                                 "nActive" = 0,
                                 "percActive" = 0
          )
          
          dataTemp$nInformed <- sum(ind_data$informed)
          dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
          
          dataTemp$nActive <- sum(ind_data$active)
          dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
          
          dataList[[t]] <- dataTemp
          
          #If all individuals are informed or an inordinately long time has elapsed, end the sim
          if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 8000) {
            break
          }
          
          #If the simulation is continuing, advance to the next time step
          t <- t + 1
        }
        
        dataCombined <- rbindlist(dataList)
        dataCombined$timeStep <- 1:nrow(dataCombined)
        
        diffusionSummary <- data.table("simID" = s, "domDist" = domDistribution,
                                       "network" = paste("dyadic", (r+1)/2, sep = "_"), "resAvail" = R,
                                       "groupAdjust" = g,
                                       "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                                       "numIndivs" = n_indivs, "numProducers" = 0,
                                       "orderDiv" = 0,
                                       "propDiv" = 0,
                                       "timeDelay" = 0,
                                       "burstiness" = 0)
        
        numProducers <- nrow(ind_data[which(ind_data$firstProd > 0),])
        orderDiv <- (1/(numProducers - 1)) * 
          sum(abs(rank(ind_data[which(ind_data$acqTime > 0  & ind_data$firstProd > 0),]$acqTime) - 
                    rank(ind_data[which(ind_data$acqTime > 0 & ind_data$firstProd > 0),]$firstProd)))
        propDiv <- sum(rank(ind_data[which(ind_data$firstProd > 0),]$acqTime) != 
                         rank(ind_data[which(ind_data$firstProd > 0),]$firstProd)) / numProducers
        timeDelay <- (1/numProducers) * sum(abs(ind_data[which(ind_data$firstProd > 0),]$acqTime - 
                                                  ind_data[which(ind_data$firstProd > 0),]$firstProd))
        
        diffusionSummary$numProducers <- numProducers
        diffusionSummary$orderDiv <- orderDiv
        diffusionSummary$propDiv <- propDiv
        diffusionSummary$timeDelay <- timeDelay
        
        intervalData <- sapply(1:(length(ind_data$acqTime) - 1), function(x) sort(ind_data$acqTime)[x+1] - sort(ind_data$acqTime)[x])
        diffusionSummary$burstiness <- (sd(intervalData)-mean(intervalData))/(sd(intervalData) + mean(intervalData))
        
        fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%s_%.02f_%s_%s.csv", run_ID, s, paste("dyad", (r+1)/2, sep = "_"), R, domDistribution, g)))
        fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%s_%.02f_%s_%s.csv", run_ID, s, paste("dyad", (r+1)/2, sep = "_"), R, domDistribution, g)))
        fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%s_%.02f_%s_%s.csv", run_ID, s, paste("dyad", (r+1)/2, sep = "_"), R, domDistribution, g)))
        
        } 
      }
      for(r in c(2,4,6)) {
        for(g in groupAdjust) {
        #Reset individual-level data to initial values
        ind_data <- ind_dataOrig
        
        domDistribution <- domDist[d]
        ind_data$distUsed <- domDistribution
        ind_data$resAvail <- R
        ind_data$network <- paste("hyper", r/2, sep = "_")
        focaldomDist <- as.vector(unlist(ind_data[domDistribution]))
        
        #Create list for holding output on each time step
        dataList <- vector("list", 8000)
        
        resources <- round(R * mean(colSums(netList[[r]])))
        
        #Set initial time
        t = 1
        
        repeat{
          
          #Identify current set of informed nodes
          informedNodes <- ind_data[which(ind_data$informed == 1),]$id
          
          #Each informed individual updates its likelihood of producing a novel behavior by:
          learningOutcomes <- hyperNetwork_diffusion(ind_data = ind_data, network = r, netList = netList, resources = resources,
                                                              informedNodes = informedNodes, domValues = focaldomDist, groupAdjustment = g)
          ind_data$active <- learningOutcomes[[1]]
          ind_data$newLearners <- learningOutcomes[[2]]
          
          for(i in 1:nrow(ind_data)) {
            if(ind_data$active[i] == 1 & ind_data$firstProd[i] == 0) {
              ind_data$firstProd[i] <- t
            }
          }
          
          #Record those individuals that acquired the trait and note the time step at which this occurred
          ind_data$acqTime <- ifelse(ind_data$newLearners == 1 & ind_data$informed == 0, t, ind_data$acqTime)
          
          #Add newly informed individuals to the list of informed individuals
          ind_data$informed <- ifelse(ind_data$informed == 0, 
                                      ifelse(ind_data$newLearners == 1, 1, 0), 1)
          
          dataTemp <- data.table("simID" = s,
                                 "domDist" = domDistribution,
                                 "network" = paste("hyper", r/2, sep = "_"),
                                 "resAvail" = R,
                                 "groupAdjust" = g,
                                 "nInformed" = 0,
                                 "percInformed" = 0,
                                 "nActive" = 0,
                                 "percActive" = 0
          )
          
          dataTemp$nInformed <- sum(ind_data$informed)
          dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
          
          dataTemp$nActive <- sum(ind_data$active)
          dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
          
          dataList[[t]] <- dataTemp
          
          #If all individuals are informed or an inordinately long time has elapsed, end the sim
          if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 8000) {
            break
          }
          
          #If the simulation is continuing, advance to the next time step
          t <- t + 1
        }
        
        dataCombined <- rbindlist(dataList)
        dataCombined$timeStep <- 1:nrow(dataCombined)
        
        diffusionSummary <- data.table("simID" = s, "domDist" = domDistribution,
                                       "network" = paste("hyper", r/2, sep = "_"), "resAvail" = R,
                                       "groupAdjust" = g,
                                       "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                                       "numIndivs" = n_indivs, "numProducers" = 0,
                                       "orderDiv" = 0,
                                       "propDiv" = 0,
                                       "timeDelay" = 0,
                                       "burstiness" = 0)
        
        numProducers <- nrow(ind_data[which(ind_data$firstProd > 0),])
        orderDiv <- (1/(numProducers - 1)) * 
          sum(abs(rank(ind_data[which(ind_data$acqTime > 0  & ind_data$firstProd > 0),]$acqTime) - 
                    rank(ind_data[which(ind_data$acqTime > 0 & ind_data$firstProd > 0),]$firstProd)))
        propDiv <- sum(rank(ind_data[which(ind_data$firstProd > 0),]$acqTime) != 
                         rank(ind_data[which(ind_data$firstProd > 0),]$firstProd)) / numProducers
        timeDelay <- (1/numProducers) * sum(abs(ind_data[which(ind_data$firstProd > 0),]$acqTime - 
                                                  ind_data[which(ind_data$firstProd > 0),]$firstProd))
        
        diffusionSummary$numProducers <- numProducers
        diffusionSummary$orderDiv <- orderDiv
        diffusionSummary$propDiv <- propDiv
        diffusionSummary$timeDelay <- timeDelay
        
        intervalData <- sapply(1:(length(ind_data$acqTime) - 1), function(x) sort(ind_data$acqTime)[x+1] - sort(ind_data$acqTime)[x])
        diffusionSummary$burstiness <- (sd(intervalData)-mean(intervalData))/(sd(intervalData) + mean(intervalData))
        
        fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%s_%.02f_%s_%s.csv", run_ID, s, paste("hyper", r/2, sep = "_"), R, domDistribution ,g)))
        fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%s_%.02f_%s_%s.csv", run_ID, s, paste("hyper", r/2, sep = "_"), R, domDistribution ,g)))
        fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%s_%.02f_%s_%s.csv", run_ID, s, paste("hyper", r/2, sep = "_"), R, domDistribution ,g)))
        
        }
      }
    }
  }
}
###############################################