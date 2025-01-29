# library(devtools)
# # install.packages("remotes")
# remotes::install_github("cran/RandomFieldsUtils")
# remotes::install_github("cran/RandomFields")
# remotes::install_github("ropensci/NLMR")

#Load packages
#source("Hypernetwork_functions.R")
source("multilayerHyperNets_simFunctions.R")
library(igraph)
#library(viridis)
library(data.table)
library(parameters)
#library(RandomFieldsUtils)
#library(RandomFields)
#library(NLMR)

bimodality <- function(skew, kurt, n) {
  beta_coef <- (skew^2+1)/(kurt + ((3 * (n-1)^2)/((n-2)*(n-3))))
  return(beta_coef)
}

setwd("C://Users/matth/Desktop/Multilayer Hypernetworks/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_indData="Sim_individual_level_data"
sim_summaryData="Sim_summary_data"
sim_subGroupSimilarity  = "Sim_subgroup_similarity_data"
sim_diffSummaryData = "Sim_diffusion_summary_data"
if(!file.exists(sim_indData)) dir.create(sim_indData)
if(!file.exists(sim_summaryData)) dir.create(sim_summaryData)
if(!file.exists(sim_diffSummaryData)) dir.create(sim_diffSummaryData)
if(!file.exists(sim_subGroupSimilarity)) dir.create(sim_subGroupSimilarity)

##Set seed for reproducibility
set.seed(01092025)

#Learning adjustment for intermediate-layer behavioral rule
#This adjusts how strongly individuals modify their production probability in response to their group mates' values
#m = 0.5

#Set social transmission rate; determines baseline probability of learning from an active neighbor in the dyadic layer
social_trans <- 0.1

c<-1.2
###################################################

nSims = 250

for(s in 1:nSims) {

##Generate individuals in family/kinship groups
ind_data <- generate_population(n_families = 25, meanFamilySize = 4, clustering = 0.075, nInitInformed = 1, clusterByFamily = FALSE)

ind_data$simID <- s

###################################################

#Generate a landscape representing relative access to an arbitrary resource (e.g., food, shelter)
# G1 <- nlm_gaussianfield(ncol = 100, nrow = 100, resolution = 1, autocorr_range = 10, mag_var = 5, 
#                         nug = 0.2, mean = 0.5, rescale = TRUE)

#Visualize resulting landscape
#plot(G1)

#Convert landscape into a matrix of resource values
#env_mat <- matrix(G1@data@values, nrow = 100, ncol = 100, byrow = TRUE)

#Assign each individual a value indicating their local resource access
#This is where having a minimum coordinate of 0.0051 is needed
#ind_data$localEnv <- sapply(1:nrow(ind_data), function(x) env_mat[round(ind_data$ind_x[x] * 100), round(ind_data$ind_y[x] * 100)])

###################################################

##Create spatial proximity network and hypergraphs
netList <- generate_latent_space_multilayer_hypergraph(ind_data = ind_data, r = c(0.15,0.2,0.25))

###############################################

#Get the mean local resources experienced by each apex-layer hyperedge
#mean_resources <- sapply(1:ncol(netList[[4]]), function(x) mean(ind_data[which(ind_data$id %in% which(netList[[4]][,x] > 0)),]$localEnv))
#ind_adjust <- sapply(1:nrow(netList[[4]]), function(x) mean(mean_resources[which(netList[[4]][x,]>0)] / ind_data$localEnv[x]))

#Identify the kin-group membership present within each upper-layer hyperedge
# familyList <- vector("list", ncol(netList[[3]]))
# for(i in 1:length(familyList)){
#   familyList[[i]] <- ind_data[which(netList[[3]][,i] > 0),]$family
# }

#This function creates adjustments based on the mean representation of relatives within an individual's hyperedges
#The 0.5 addition means that if half of the hyperedge (group) membership are kin, then there is no adjustment (i.e., adjustment = 1)
#Otherwise, if kin membership is over half, the adjustment is > 1; less than half means its < 1
#ind_adjust2 <- sapply(1:nrow(netList[[3]]), function(y) mean(sapply(which(netList[[3]][y,]>0), function(x) (sum(familyList[[x]]==ind_data[y,]$family) / 
#                                                                                                              length(familyList[[x]])) + 0.5)))
#ind_data$resourceAdjust <- ind_adjust
#ind_data$familyAdjust <- ind_adjust2

ind_dataOrig <- ind_data

#ruleSet <- c("deltaSmartest", "deltaCombined", "dampMin", "deltaAll")
netType <- c("dyadic", "higherOrder")

for(r in 1:length(netType)) {
  
  #Reset individual-level data to initial values
  ind_data <- ind_dataOrig
  
  #adjustRule <- ruleSet[r]
  
  ind_data$rule <- adjustRule <- "domRank"
  
  ind_data$netType <- net <- netType[r]
  
  ind_data$firstProd <- 0
  
#Create list for holding output on each time step
dataList <- vector("list", 5000)

#Set initial time
t = 1

repeat{
  
  #Identify current set of informed nodes
  informedNodes <- ind_data[which(ind_data$informed == 1),]$id
  
  #Each informed individual updates its likelihood of producing a novel behavior by: 
  #first, its access to resources relative to others in the apex-layer hyperedges to which it belongs, and
  #second, the relative proportion of kin within its upper-layer hyperedges
  #For now I cap out an individual's likelhood of production at 1 if the resulting adjustment would put it above 1
  # for(i in informedNodes) {
  #   ind_data$produce[i] <- min(1, ind_adjust[i] * ind_data$produce[i])
  #   ind_data$produce[i] <- min(1, ind_adjust2[i] *  ind_data$produce[i])
  # }
  
  #Keeping this note here for now, but the code above addresses points 2 and 3; point 1 may or may not be an issue
  #Three issues here that I can see: (1) using mean values of resources for each group produces quite a narrow range, with most groups hovering around 0.5
  #ranging from 0.4 to ~0.65; (2) modifying individual's produce likelihoods here can raise them above 1; perhaps need to re-normalize?; and (3) all individuals
  #adjust their values here, but may only want to do this for individuals that are informed.
  
  produceTemp <- produceByDomRank(ind_data = ind_data, netList = netList, domSteepness = c, netType = net)
  
  #Can implement one of deltaSmartest, deltaAll, dampMin, deltaCombined
  #produceTemp <- adjust_production_probability(ind_data = ind_data, netList = netList, rule = adjustRule)
  
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
  ind_data$acqProb <- sapply(1:nrow(ind_data), function(x) social_trans * sum(netList[[1]][,x] * ind_data$informed) / sum(netList[[1]][,x]))
  ind_data$acqProb <- ifelse(is.na(ind_data$acqProb), 0, ind_data$acqProb)
  
  #Record those individuals that acquired the trait and note the time step at which this occurred
  acqVect <- runif(nrow(ind_data), min = 0, max = 1)
  ind_data$learned <- 0
  ind_data$learned <- ifelse(ind_data$acqProb > acqVect, 1, 0)
  ind_data$acqTime <- ifelse(ind_data$learned == 1 & ind_data$informed == 0, t, ind_data$acqTime)
  
  #Add newly informed individuals to the list of informed individuals
  ind_data$informed <- ifelse(ind_data$informed == 0, 
                              ifelse(ind_data$learned == 1, 1, 0), 1)
  
  n_indivs <- nrow(ind_data)
  
  dataTemp <- data.table("simID" = s, 
                         "rule" = adjustRule,
                         "nInformed" = 0,
                         "percInformed" = 0,
                         "nActive" = 0,
                         "percActive" = 0,
                         "meanProd_I" = 0,
                         "varProd_I" = 0,
                         "skewProd_I" = 0,
                         "kurtProd_I" = 0,
                         "bimodProd_I" = 0
  )
  
  dataTemp$nInformed <- sum(ind_data$informed)
  dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
  
  dataTemp$nActive <- sum(ind_data$active)
  dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
  
  dataTemp$meanProd_I <- mean(ind_data[which(ind_data$informed == 1),]$produce)
  dataTemp$varProd_I <- var(ind_data[which(ind_data$informed == 1),]$produce)
  if(dataTemp$nInformed >= 3) {
    dataTemp$skewProd_I <- skewness(ind_data[which(ind_data$informed == 1),]$produce)[[1]]
  } else {
    dataTemp$skewProd_I <- NA
  } 
  if(dataTemp$nInformed >= 4) {
    dataTemp$kurtProd_I <- kurtosis(ind_data[which(ind_data$informed == 1),]$produce)[[1]]
    dataTemp$bimodProd_I <- bimodality(skew = dataTemp$skewProd_I, kurt = dataTemp$kurtProd_I, n = n_indivs)
  } else{
    dataTemp$kurtProd_I <- NA
    dataTemp$bimodProd_I <- NA
  }
  
  dataList[[t]] <- dataTemp
  
  #If all individuals are informed or an inordinately long time has elapsed, end the sim
  if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 5000) {
    #sum(ind_data$informed) == nrow(ind_data) | t == 5000) {
    break
  }
  
  #If the simulation is continuing, advance to the next time step
  t <- t + 1
}
}

dataCombined <- rbindlist(dataList)
dataCombined$timeStep <- 1:nrow(dataCombined)
dataCombined$netType <- net

diffusionSummary <- data.table("simID" = s, "rule" = adjustRule, "netType" = net,
                               "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                               "orderDiv" = (1/(n_indivs - 1))*sum(abs(rank(ind_data[which(ind_data$acqTime > 0),]$acqTime) - rank(ind_data[which(ind_data$acqTime > 0),]$firstProd))),
                               "propDiv" = sum(rank(ind_data$acqTime) != rank(ind_data$firstProd))/n_indivs,
                               "timeDelay" = (1/n_indivs) * sum(abs(ind_data$acqTime - ind_data$firstProd)),
                               "burstiness" = 0)

intervalData <- sapply(1:(length(ind_data$acqTime) - 1), function(x) sort(ind_data$acqTime)[x+1] - sort(ind_data$acqTime)[x])
diffusionSummary$burstiness <- (sd(intervalData)-mean(intervalData))/(sd(intervalData) + mean(intervalData))

subGroupData <- data.frame("simID" = s, "rule" = adjustRule, 
                           "SubID" = 1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), 
                           "similarity" = sapply(1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), function(x) mean(as.numeric(dist(ind_data[which(ind_data$id %in% as.vector(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))[[x]])),]$produce)))),
                           "subGroupSize" = sapply(1:length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))), function(x) length(max_cliques(graph_from_adjacency_matrix(netList[[1]], mode = "max"))[[x]])))

fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%s.csv", run_ID, s, net)))
fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%s.csv", run_ID, s, net)))
fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%s.csv", run_ID, s, net)))
fwrite(subGroupData, file = file.path(sim_subGroupSimilarity, sprintf("simData_%s_%.01i_%s.csv", run_ID, s, net)))

}

# 
# plot(dataCombined$timeStep, dataCombined$percInformed)
# plot(dataCombined$timeStep, dataCombined$meanProd)
# plot(dataCombined$timeStep, dataCombined$varProd)
# plot(dataCombined$timeStep, dataCombined$skewProd)
# plot(dataCombined$timeStep, dataCombined$kurtProd)
# plot(dataCombined$timeStep, dataCombined$bimodProd)
# plot(dataCombined$timeStep, dataCombined$meanProd_I)
# plot(dataCombined$timeStep, dataCombined$varProd_I)
# plot(dataCombined$timeStep, dataCombined$skewProd_I)
# plot(dataCombined$timeStep, dataCombined$kurtProd_I)
# plot(dataCombined$timeStep, dataCombined$bimodProd_I)
# plot(dataCombined$timeStep, dataCombined$resCor)
# plot(dataCombined$timeStep, dataCombined$resCor_I)
# plot(dataCombined$timeStep, dataCombined$famCor)
# plot(dataCombined$timeStep, dataCombined$famCor_I)
# plot(familyData$familySize, familyData$similarity)
# hist(familyData$similarity, breaks = 20)
# hist(subGroupData$similarity, breaks = 20)
# plot(ind_data$resourceAdjust, ind_data$produce)
# plot(ind_data$familyAdjust, ind_data$produce)

###############################################

