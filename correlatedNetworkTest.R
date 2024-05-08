source("Hypernetwork_functions.R")
source("dataImportAndAggregation.R")
library(gdata)
library(data.table)
library(doParallel)
library(igraph)
library(dplyr)
library(purrr)

registerDoParallel(cores = 20)

# generate_correlated_hypergraph <- function(hypergraph, rewireP) {
#   numEdges <- ncol(hypergraph)
#   rewire <- runif(numEdges, min = 0, max = 1)
#   rewire <- ifelse(rewire <= rewireP, 1, 0)
#   rewireList <- which(rewire == 1)
#   if(length(rewireList) > 0) {
#     rewireList <- resample(rewireList, length(rewireList))
#     for(i in rewireList) {
#       swapEdge <- sample(which(!(1:numEdges %in% i)),1)
#       Hi <- resample(which(hypergraph[,i] > 0),1)
#       Hx <- resample(which(hypergraph[,swapEdge] > 0),1)
#       hypergraph[Hi,i] <- 0
#       hypergraph[Hx,i] <- 1
#     }
#   }
#   return(hypergraph)
# }
# 
# setwd("C://Users/matth/Desktop/")
# 
# #Create folder in which to store simulation results
# run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
# sim_details="hyperGraph_Sequences_lowCorr_May24"
# if(!file.exists(sim_details)) dir.create(sim_details)
# 
# set.seed(05082024)
# 
# for(s in 1:100){
#   
#   H1 <- generate_ER_hypergraph(v = 500, m = 500, p = 0.01)
#   
#   for(t in 1:20) {
#     
#     fwrite(H1, file = file.path(sim_details, sprintf("simData_%s_%.01i_%.01i.csv", run_ID, s, t)))
#     H1 <- generate_correlated_hypergraph(hypergraph = H1, rewireP = 0.5)
#     
#   }
# }

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details_low="saTBCData_lowCorr_May24"
sim_details_high="saTBCData_highCorr_May24"
if(!file.exists(sim_details_low)) dir.create(sim_details_low)
if(!file.exists(sim_details_high)) dir.create(sim_details_high)

highCorList <- importCSVs(path = "~/scratch/Temporal_HyperNets/Run1/hyperGraph_Sequences_highCorr_May24/")
lowCorList <- importCSVs(path = "~/scratch/Temporal_HyperNets/Run1/hyperGraph_Sequences_lowCorr_May24/")

highCorIncidMats <- vector("list", 100)

for(s in 1:length(highCorList)) {
i <- ceiling(s/20)
t <- ifelse((s %% 20) == 0, 20, s %% 20)

highCorIncidMats[[i]][[t]] <- highCorList[[s]]

}

lowCorIncidMats <- vector("list", 100)

for(s in 1:length(lowCorList)) {
  i <- ceiling(s/20)
  t <- ifelse((s %% 20) == 0, 20, s %% 20)
  
  lowCorIncidMats[[i]][[t]] <- lowCorList[[s]]
  
}

alphaValues <- c(0.01, 0.5, 0.99)
shortWindowVect <- c(1, 3, 5)
longWindowVect <- c(1, 10, 20)

foreach(s = 1:length(lowCorIncidMats)) %dopar% {
  shortWindowStack <- vector("list", length(shortWindowVect))
  longWindowStack <- vector("list", length(longWindowVect))
  
  for(t in 1:length(shortWindowVect)) {
    shortWindowStack[[t]] <- lowCorIncidMats[[s]][[shortWindowVect[t]]]
    longWindowStack[[t]] <- lowCorIncidMats[[s]][[longWindowVect[t]]]
  }
  #nCores <- detectCores() - 1
  for(a in alphaValues) {
    saTBC <-get_sa_temporalBC(hypergraphList = shortWindowStack, timeStamps = shortWindowVect, smax = 4, windowLength = 20, focalTimeStamp = 20, alpha = a, 
                      normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s, 
                           "rewireP" = 0.5, 
                           "window" = 3,
                           "alpha" = a, 
                           "ID" = names(saTBC[[2]]), 
                           "saTBC" = saTBC[[2]])
    write.csv(dataTemp, file = file.path(sim_details_low, sprintf("simData_%s_%.01i_%.02f_%.01i_%.02f.csv", run_ID, s, 0.5, 3, a)))
    
    saTBC2 <-get_sa_temporalBC(hypergraphList = longWindowStack, timeStamps = longWindowVect, smax = 4, windowLength = 20, focalTimeStamp = 20, alpha = a, 
                              normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s, 
                           "rewireP" = 0.5, 
                           "window" = 10,
                           "alpha" = a, 
                           "ID" = names(saTBC2[[2]]), 
                           "saTBC" = saTBC2[[2]])
    write.csv(dataTemp, file = file.path(sim_details_low, sprintf("simData_%s_%.01i_%.02f_%.02i_%.02f.csv", run_ID, s, 0.5, 10, a)))
  }
  
  shortWindowStack <- vector("list", length(shortWindowVect))
  longWindowStack <- vector("list", length(longWindowVect))
  
  for(t in 1:length(shortWindowVect)) {
    shortWindowStack[[t]] <- highCorIncidMats[[s]][[shortWindowVect[t]]]
    longWindowStack[[t]] <- highCorIncidMats[[s]][[longWindowVect[t]]]
  }
  #nCores <- detectCores() - 1
  for(a in alphaValues) {
    saTBC <-get_sa_temporalBC(hypergraphList = shortWindowStack, timeStamps = shortWindowVect, smax = 4, windowLength = 20, focalTimeStamp = 20, alpha = a, 
                              normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s, 
                           "rewireP" = 0.1, 
                           "window" = 3,
                           "alpha" = a, 
                           "ID" = names(saTBC[[2]]), 
                           "saTBC" = saTBC[[2]])
    write.csv(dataTemp, file = file.path(sim_details_high, sprintf("simData_%s_%.01i_%.02f_%.01i_%.02f.csv", run_ID, s, 0.1, 3, a)))
    
    saTBC2 <-get_sa_temporalBC(hypergraphList = longWindowStack, timeStamps = longWindowVect, smax = 4, windowLength = 20, focalTimeStamp = 20, alpha = a, 
                               normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s, 
                           "rewireP" = 0.1, 
                           "window" = 10,
                           "alpha" = a, 
                           "ID" = names(saTBC2[[2]]), 
                           "saTBC" = saTBC2[[2]])
    write.csv(dataTemp, file = file.path(sim_details_high, sprintf("simData_%s_%.01i_%.02f_%.02i_%.02f.csv", run_ID, s, 0.1, 10, a)))
  }
  
}


