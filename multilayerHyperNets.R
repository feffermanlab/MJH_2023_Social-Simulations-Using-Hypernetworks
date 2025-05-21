
#Keeping here for now in case I want to include resource distributions down the road
# library(devtools)
# # install.packages("remotes")
# remotes::install_github("cran/RandomFieldsUtils")
# remotes::install_github("cran/RandomFields")
# remotes::install_github("ropensci/NLMR")
#library(RandomFieldsUtils)
#library(RandomFields)
#library(NLMR)

#Load packages
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/multilayerHyperNets_simFunctions.R", chdir = TRUE)
#source("multilayerHyperNets_simFunctions.R")
library(igraph)
library(data.table)
library(foreach)
library(doParallel)
library(stringr)

registerDoParallel(cores = 20)

# bimodality <- function(skew, kurt, n) {
#   beta_coef <- (skew^2+1)/(kurt + ((3 * (n-1)^2)/((n-2)*(n-3))))
#   return(beta_coef)
# }

#setwd("C://Users/matth/Desktop/Multilayer Hypernetworks/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_indData="Sim_individual_level_data"
sim_summaryData="Sim_summary_data"
sim_subGroupSimilarity  = "Sim_subgroup_similarity_data"
sim_diffSummaryData = "Sim_diffusion_summary_data"
sim_netData_dyadic = "Sim_dyadic_network"
sim_netData_HE1 = "Sim_higherOrder_network_1"
sim_netData_HE2 = "Sim_higherOrder_network_2"
sim_netData_HE3 = "Sim_higherOrder_network_3"

# sim_indData="Sim_individual_level_data_Inc"
# sim_summaryData="Sim_summary_data_Inc"
# sim_subGroupSimilarity  = "Sim_subgroup_similarity_data_Inc"
# sim_diffSummaryData = "Sim_diffusion_summary_data_Inc"
# sim_netData_dyadic = "Sim_dyadic_network_Inc"
# sim_netData_HE1 = "Sim_higherOrder_network_1_Inc"
# sim_netData_HE2 = "Sim_higherOrder_network_2_Inc"
# sim_netData_HE3 = "Sim_higherOrder_network_3_Inc"

if(!file.exists(sim_indData)) dir.create(sim_indData)
if(!file.exists(sim_summaryData)) dir.create(sim_summaryData)
if(!file.exists(sim_diffSummaryData)) dir.create(sim_diffSummaryData)
if(!file.exists(sim_subGroupSimilarity)) dir.create(sim_subGroupSimilarity)
if(!file.exists(sim_netData_dyadic)) dir.create(sim_netData_dyadic)
if(!file.exists(sim_netData_HE1)) dir.create(sim_netData_HE1)
if(!file.exists(sim_netData_HE2)) dir.create(sim_netData_HE2)
if(!file.exists(sim_netData_HE3)) dir.create(sim_netData_HE3)

##Set seed for reproducibility
#set.seed(01092025)
#set.seed(03042025)
#set.seed(04102025)
set.seed(5212025)

#Set social transmission rate; determines baseline probability of learning from an active neighbor in the dyadic layer
social_trans <- 0.1

#Possible radii for determinig subgroups
rad <- c(0.15, 0.25, 0.35)

#Possible functions for determining resopnse to dominant individuals
domResp <- c(domResponse_Linear_HE, domResponse_Sigmoid2_HE, domResponse_Linear_dyadic, domResponse_Sigmoid2_dyadic,
             domResponse_wLinear_HE, domResponse_wSigmoid2_HE, domResponse_wLinear_dyadic, domResponse_wSigmoid2_dyadic)

# domResp <- c(domResponse_Linear_HE_Increasing, domResponse_Sigmoid2_HE_Increasing, domResponse_Linear_dyadic_Increasing, domResponse_Sigmoid2_dyadic_Increasing, 
#              domResponse_wLinear_HE_Increasing, domResponse_wSigmoid2_HE_Increasing, domResponse_wLinear_dyadic_Increasing, domResponse_wSigmoid2_dyadic_Increasing)

domRespNames <- c("Linear_HE", "Sigmoid2_HE", "Linear_dyadic", "Sigmoid2_dyadic",
                  "wLinear_HE", "wSigmoid2_HE", "wLinear_dyadic", "wSigmoid2_dyadic")

#Distributions for dominance ranks
domDist <- c("domUni", "domExp")

###################################################

#Set number of sims per condition
nSims = 250

foreach(s = 1:nSims) %dopar% {
#for(s in 1:nSims) {
  
  ##Generate individuals in family/kinship groups
  #For the moment, not using families in model, but keeping in case of expanding
  ind_data <- generate_population(n_families = 25, meanFamilySize = 4, clustering = 0.075, 
                                  nInitInformed = 1, clusterByFamily = FALSE, minRadius = rad[1])
  n_indivs <- nrow(ind_data)
  ind_data$simID <- s
  
  ind_data$domUni <- runif(nrow(ind_data), min = 0, max = 1)
  domExp <- sort(rexp(nrow(ind_data), rate = 2))
  domExp <- sapply(rank(ind_data$domUni), function(x) domExp[x])
  ind_data$domExp <- (domExp - min(domExp))/(max(domExp) - min(domExp))
  ind_data$domExp <- ifelse(ind_data$domExp == 0, min(ind_data[which(ind_data$domExp > 0),]$domExp - 0.0001), ind_data$domExp)
  
  ind_dataOrig <- ind_data
  
  ##Create spatial proximity network and hypergraphs
  netList <- generate_latent_space_multilayer_hypergraph(ind_data = ind_data, r = rad)
  
  write.csv(netList[[1]], file.path(sim_netData_dyadic, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[2]], file.path(sim_netData_HE1, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[3]], file.path(sim_netData_HE2, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[4]], file.path(sim_netData_HE3, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  
  for(d in 1:length(domDist)) {
    
    for(p in 1:length(domResp)) {
      
      if(str_split(domRespNames[p], "_")[[1]][2] == "dyadic") {
        
      #Reset individual-level data to initial values
      ind_data <- ind_dataOrig
      
      domDistribution <- domDist[d]
      
      ind_data$domResponse <- domRespNames[p]
      ind_data$distUsed <- domDistribution
      ind_data$groupRadius <- 0
      
      focaldomDist <- as.vector(unlist(ind_data[domDistribution]))
      domResponseFunction <- domResp[[p]]
      radiusForDomFunction <- 1
      
      ind_data$firstProd <- 0
      
      #Create list for holding output on each time step
      dataList <- vector("list", 20000)
      
      #Set initial time
      t = 1
      
      repeat{
        
        #Identify current set of informed nodes
        informedNodes <- ind_data[which(ind_data$informed == 1),]$id
        
        #Each informed individual updates its likelihood of producing a novel behavior by:
        produceTemp <- domResponseFunction(ind_data = ind_data, netList = netList, radius = radiusForDomFunction, informedNodes = informedNodes, domValues = focaldomDist)
        ind_data$produce <- produceTemp
        
        #Determine which informed nodes are producing the novel behavior this time step
        probVect <- runif(nrow(ind_data), min = 0, max = 1)
        ind_data$active <- 0
        
        for(i in informedNodes) {
          ind_data$active[i] <- ifelse(ind_data$produce[i] >  probVect[i], 1, 0)
        }
        
        for(i in 1:nrow(ind_data)) {
          if(ind_data$active[i] == 1 & ind_data$firstProd[i] == 0) {
            ind_data$firstProd[i] <- t
          }
        }
        
        #Determine individuals' likelihood of learning the trait based on their relative connection strength to informed individuals in the dyadic layer
        ind_data$acqProb <- sapply(1:nrow(ind_data), function(x) social_trans * sum(netList[[1]][,x] * ind_data$active) / sum(netList[[1]][,x]))
        ind_data$acqProb <- ifelse(is.na(ind_data$acqProb), 0, ind_data$acqProb)
        
        #Record those individuals that acquired the trait and note the time step at which this occurred
        acqVect <- runif(nrow(ind_data), min = 0, max = 1)
        ind_data$learned <- 0
        ind_data$learned <- ifelse(ind_data$acqProb > acqVect, 1, 0)
        ind_data$acqTime <- ifelse(ind_data$learned == 1 & ind_data$informed == 0, t, ind_data$acqTime)
        
        #Add newly informed individuals to the list of informed individuals
        ind_data$informed <- ifelse(ind_data$informed == 0, 
                                    ifelse(ind_data$learned == 1, 1, 0), 1)
        
        dataTemp <- data.table("simID" = s,
                               "domDist" = domDistribution,
                               "domResponse" = domRespNames[p],
                               "groupRadius" = rad[r],
                               "nInformed" = 0,
                               "percInformed" = 0,
                               "nActive" = 0,
                               "percActive" = 0,
                               "meanProd_I" = 0,
                               "varProd_I" = 0
        )
        
        dataTemp$nInformed <- sum(ind_data$informed)
        dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
        
        dataTemp$nActive <- sum(ind_data$active)
        dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
        
        dataTemp$meanProd_I <- mean(ind_data[which(ind_data$informed == 1),]$produce)
        dataTemp$varProd_I <- var(ind_data[which(ind_data$informed == 1),]$produce)
        
        dataList[[t]] <- dataTemp
        
        #If all individuals are informed or an inordinately long time has elapsed, end the sim
        if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 20000) {
          break
        }
        
        #If the simulation is continuing, advance to the next time step
        t <- t + 1
      }
      
      dataCombined <- rbindlist(dataList)
      dataCombined$timeStep <- 1:nrow(dataCombined)
      
      diffusionSummary <- data.table("simID" = s, "domDist" = domDistribution,
                                     "domResponse" = domRespNames[p], "groupRadius" = rad[r],
                                     "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                                     "numIndivs" = n_indivs, "numProducers" = 0, "numInformed" = 0,
                                     "orderDiv" = 0,
                                     "propDiv" = 0,
                                     "timeDelay" = 0,
                                     "burstiness" = 0)
      
      numProducers <- nrow(ind_data[which(ind_data$firstProd > 0),])
      numInformed <- nrow(ind_data[which(ind_data$informed > 0),])
      orderDiv <- (1/(numProducers - 1)) * 
        sum(abs(rank(ind_data[which(ind_data$acqTime > 0  & ind_data$firstProd > 0),]$acqTime) - 
                  rank(ind_data[which(ind_data$acqTime > 0 & ind_data$firstProd > 0),]$firstProd)))
      propDiv <- sum(rank(ind_data[which(ind_data$firstProd > 0),]$acqTime) != 
                       rank(ind_data[which(ind_data$firstProd > 0),]$firstProd)) / numProducers
      timeDelay <- (1/numProducers) * sum(abs(ind_data[which(ind_data$firstProd > 0),]$acqTime - 
                                                ind_data[which(ind_data$firstProd > 0),]$firstProd))
      
      diffusionSummary$numInformed <- numInformed
      diffusionSummary$numProducers <- numProducers
      diffusionSummary$orderDiv <- orderDiv
      diffusionSummary$propDiv <- propDiv
      diffusionSummary$timeDelay <- timeDelay
      
      intervalData <- sapply(1:(length(ind_data$acqTime) - 1), function(x) sort(ind_data$acqTime)[x+1] - sort(ind_data$acqTime)[x])
      diffusionSummary$burstiness <- (sd(intervalData)-mean(intervalData))/(sd(intervalData) + mean(intervalData))
      
      subGroupData <- data.frame("simID" = s, "domDist" = domDistribution,
                                 "domResponse" = domRespNames[p], "groupRadius" = rad[r],
                                 "SubID" = 1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), 
                                 "similarity" = sapply(1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), function(x) mean(as.numeric(dist(ind_data[which(ind_data$id %in% as.vector(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))[[x]])),]$produce)))),
                                 "subGroupSize" = sapply(1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), function(x) length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))[[x]])))
      
      fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
      fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
      fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
      fwrite(subGroupData, file = file.path(sim_subGroupSimilarity, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
      
      } else{
      if(str_split(domRespNames[p], "_")[[1]][2] == "HE") {
      for(r in 1:length(rad)) {

        #Reset individual-level data to initial values
        ind_data <- ind_dataOrig
        
        domDistribution <- domDist[d]
        
        ind_data$domResponse <- domRespNames[p]
        ind_data$distUsed <- domDistribution
        ind_data$groupRadius <- rad[r]
        
        focaldomDist <- as.vector(unlist(ind_data[domDistribution]))
        domResponseFunction <- domResp[[p]]
        radiusForDomFunction <- r + 1
        
        ind_data$firstProd <- 0
        
        #Create list for holding output on each time step
        dataList <- vector("list", 20000)
        
        #Set initial time
        t = 1
        
        repeat{
          
          #Identify current set of informed nodes
          informedNodes <- ind_data[which(ind_data$informed == 1),]$id
          
          #Each informed individual updates its likelihood of producing a novel behavior by:
          produceTemp <- domResponseFunction(ind_data = ind_data, netList = netList, radius = radiusForDomFunction, informedNodes = informedNodes, domValues = focaldomDist)
          ind_data$produce <- produceTemp
          
          #Determine which informed nodes are producing the novel behavior this time step
          probVect <- runif(nrow(ind_data), min = 0, max = 1)
          ind_data$active <- 0
          
          for(i in informedNodes) {
            ind_data$active[i] <- ifelse(ind_data$produce[i] >  probVect[i], 1, 0)
          }
          
          for(i in 1:nrow(ind_data)) {
            if(ind_data$active[i] == 1 & ind_data$firstProd[i] == 0) {
              ind_data$firstProd[i] <- t
            }
          }
          
          #Determine individuals' likelihood of learning the trait based on their relative connection strength to informed individuals in the dyadic layer
          ind_data$acqProb <- sapply(1:nrow(ind_data), function(x) social_trans * sum(netList[[1]][,x] * ind_data$active) / sum(netList[[1]][,x]))
          ind_data$acqProb <- ifelse(is.na(ind_data$acqProb), 0, ind_data$acqProb)
          
          #Record those individuals that acquired the trait and note the time step at which this occurred
          acqVect <- runif(nrow(ind_data), min = 0, max = 1)
          ind_data$learned <- 0
          ind_data$learned <- ifelse(ind_data$acqProb > acqVect, 1, 0)
          ind_data$acqTime <- ifelse(ind_data$learned == 1 & ind_data$informed == 0, t, ind_data$acqTime)
          
          #Add newly informed individuals to the list of informed individuals
          ind_data$informed <- ifelse(ind_data$informed == 0, 
                                      ifelse(ind_data$learned == 1, 1, 0), 1)
          
          dataTemp <- data.table("simID" = s,
                                 "domDist" = domDistribution,
                                 "domResponse" = domRespNames[p],
                                 "groupRadius" = rad[r],
                                 "nInformed" = 0,
                                 "percInformed" = 0,
                                 "nActive" = 0,
                                 "percActive" = 0,
                                 "meanProd_I" = 0,
                                 "varProd_I" = 0
          )
          
          dataTemp$nInformed <- sum(ind_data$informed)
          dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
          
          dataTemp$nActive <- sum(ind_data$active)
          dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
          
          dataTemp$meanProd_I <- mean(ind_data[which(ind_data$informed == 1),]$produce)
          dataTemp$varProd_I <- var(ind_data[which(ind_data$informed == 1),]$produce)
          
          dataList[[t]] <- dataTemp
          
          #If all individuals are informed or an inordinately long time has elapsed, end the sim
          if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 20000) {
            break
          }
          
          #If the simulation is continuing, advance to the next time step
          t <- t + 1
        }
        
        dataCombined <- rbindlist(dataList)
        dataCombined$timeStep <- 1:nrow(dataCombined)
        
        diffusionSummary <- data.table("simID" = s, "domDist" = domDistribution,
                                       "domResponse" = domRespNames[p], "groupRadius" = rad[r],
                                       "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                                       "numIndivs" = n_indivs, "numProducers" = 0, "numInformed" = 0,
                                       "orderDiv" = 0,
                                       "propDiv" = 0,
                                       "timeDelay" = 0,
                                       "burstiness" = 0)
        
        numInformed <- nrow(ind_data[which(ind_data$informed > 0),])
        numProducers <- nrow(ind_data[which(ind_data$firstProd > 0),])
        orderDiv <- (1/(numProducers - 1)) * 
          sum(abs(rank(ind_data[which(ind_data$acqTime > 0  & ind_data$firstProd > 0),]$acqTime) - 
                    rank(ind_data[which(ind_data$acqTime > 0 & ind_data$firstProd > 0),]$firstProd)))
        propDiv <- sum(rank(ind_data[which(ind_data$firstProd > 0),]$acqTime) != 
                         rank(ind_data[which(ind_data$firstProd > 0),]$firstProd)) / numProducers
        timeDelay <- (1/numProducers) * sum(abs(ind_data[which(ind_data$firstProd > 0),]$acqTime - 
                                                  ind_data[which(ind_data$firstProd > 0),]$firstProd))
        
        diffusionSummary$numInformed <- numInformed
        diffusionSummary$numProducers <- numProducers
        diffusionSummary$orderDiv <- orderDiv
        diffusionSummary$propDiv <- propDiv
        diffusionSummary$timeDelay <- timeDelay
        
        intervalData <- sapply(1:(length(ind_data$acqTime) - 1), function(x) sort(ind_data$acqTime)[x+1] - sort(ind_data$acqTime)[x])
        diffusionSummary$burstiness <- (sd(intervalData)-mean(intervalData))/(sd(intervalData) + mean(intervalData))
        
        subGroupData <- data.frame("simID" = s, "domDist" = domDistribution,
                                   "domResponse" = domRespNames[p], "groupRadius" = rad[r],
                                   "SubID" = 1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), 
                                   "similarity" = sapply(1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), function(x) mean(as.numeric(dist(ind_data[which(ind_data$id %in% as.vector(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))[[x]])),]$produce)))),
                                   "subGroupSize" = sapply(1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), function(x) length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))[[x]])))
        
        fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        fwrite(subGroupData, file = file.path(sim_subGroupSimilarity, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        
      }
      }
      }
    }
  }
}

###############################################

