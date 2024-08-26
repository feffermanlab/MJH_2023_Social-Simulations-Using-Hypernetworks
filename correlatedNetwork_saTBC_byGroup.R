source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
library(data.table)
library(doParallel)
library(igraph)
library(dplyr)
library(purrr)
library(stringr)

registerDoParallel(cores = 20)

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details_low="saTBCData_byGroup_lowCorr"
sim_details_high="saTBCData_byGroup_highCorr"
if(!file.exists(sim_details_low)) dir.create(sim_details_low)
if(!file.exists(sim_details_high)) dir.create(sim_details_high)

highCorList <- importCSVs(path = "~/scratch/Temporal_HyperNets/Run2/hyperGraph_Sequences_highCorr_Aug23/")
lowCorList <- importCSVs(path = "~/scratch/Temporal_HyperNets/Run2/hyperGraph_Sequences_lowCorr_Aug23/")

highCorIncidMats <- vector("list", 100)

for(s in 1:length(highCorList)) {
i <- ceiling(s/20)
t <- ifelse((s %% 20) == 0, 20, s %% 20)

highCorIncidMats[[i]][[t]] <- t(highCorList[[s]])

}

lowCorIncidMats <- vector("list", 100)

for(s in 1:length(lowCorList)) {
  i <- ceiling(s/20)
  t <- ifelse((s %% 20) == 0, 20, s %% 20)
  lowCorIncidMats[[i]][[t]] <- t(lowCorList[[s]])
}

alphaValues <- c(0.01, 0.99)
shortWindowLength = 3
longWindowLength = 10
focalTimePoint = 20
shortWindowVect = seq(from = focalTimePoint - (shortWindowLength - 1), to = focalTimePoint)
longWindowVect = seq(from = focalTimePoint - (longWindowLength - 1), to = focalTimePoint)

foreach(s = 1:length(lowCorIncidMats)) %dopar% {

  for(a in alphaValues) {
    saTBC <-get_sa_temporalBC(hypergraphList = lowCorIncidMats[[s]], timeStamps = shortWindowVect, smax = 12, 
                              windowLength = 20, focalTimeStamp = 20, alpha = a, normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s, 
                           "rewireP" = 0.9, 
                           "window" = 3,
                           "alpha" = a, 
                           "ID" = names(saTBC[[2]])) 
    write.csv(cbind(dataTemp,saTBC[[1]]), file = file.path(sim_details_low, sprintf("simData_%s_%.01i_%.02f_%.01i_%.02f.csv", run_ID, s, 0.5, 3, a)))
    
    saTBC2 <-get_sa_temporalBC(hypergraphList = lowCorIncidMats[[s]], timeStamps = longWindowVect, smax = 12,
                               windowLength = 20, focalTimeStamp = 20, alpha = a, normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s, 
                           "rewireP" = 0.9, 
                           "window" = 10,
                           "alpha" = a, 
                           "ID" = names(saTBC2[[2]]))
    write.csv(cbind(dataTemp,saTBC2[[1]]), file = file.path(sim_details_low, sprintf("simData_%s_%.01i_%.02f_%.02i_%.02f.csv", run_ID, s, 0.5, 10, a)))
  }

  for(a in alphaValues) {
    saTBC <-get_sa_temporalBC(hypergraphList = highCorIncidMats[[s]], timeStamps = shortWindowVect, smax = 12, 
                              windowLength = 20, focalTimeStamp = 20, alpha = a, normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s,
                           "rewireP" = 0.1,
                           "window" = 3,
                           "alpha" = a,
                           "ID" = names(saTBC[[2]]))
    write.csv(cbind(dataTemp, saTBC[[1]]), file = file.path(sim_details_high, sprintf("simData_%s_%.01i_%.02f_%.01i_%.02f.csv", run_ID, s, 0.1, 3, a)))

    saTBC2 <-get_sa_temporalBC(hypergraphList = highCorIncidMats[[s]], timeStamps = longWindowVect, smax = 12, 
                               windowLength = 20, focalTimeStamp = 20, alpha = a, normalized = TRUE, method = "serial", nCores = NULL)
    dataTemp <- data.table("simID" = s,
                           "rewireP" = 0.1,
                           "window" = 10,
                           "alpha" = a,
                           "ID" = names(saTBC2[[2]]))
    write.csv(cbind(dataTemp,saTBC2[[1]]), file = file.path(sim_details_high, sprintf("simData_%s_%.01i_%.02f_%.02i_%.02f.csv", run_ID, s, 0.1, 10, a)))
  }
}
