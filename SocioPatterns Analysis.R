
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)

library(data.table)
library(purrr)
library(igraph)
library(dplyr)
library(foreach)
library(doParallel)
library(stringr)

# jaccard <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }

#Create folders in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sTBC_results = "sTBC_SFHH_results"
if(!file.exists(edge_lists)) dir.create(sTBC_results)

#Import incidence matrices

incidMatList <- importCSVs(path = "/home/mhasenja/scratch/Temporal_HyperNets/SocioPatterns/Incidence Matrices/SFHH/")

for(i in 1:length(incidMatList)) {
  rownames(incidMatList[[i]]) <- incidMatList[[i]][,1]
  incidMatList[[i]] <- incidMatList[[i]][,-1]
  incidMatList[[i]] <- t(incidMatList[[i]])
}
timeStamps <- c(1:2)

# indiv_saTBC_alpha0.5 <- get_sa_temporalBC(hypergraphList = incidMatList, timeStamps = timeStamps, smax = 25, windowLength = 2, 
#                                           focalTimeStamp = 2, alpha = 0.5, normalized = TRUE, method = "parallel", nCores = 11)
# 
# indiv_saTBC_alpha0.99 <- get_sa_temporalBC(hypergraphList = incidMatList, timeStamps = timeStamps, smax = 25, windowLength = 2, 
#                                            focalTimeStamp = 2, alpha = 0.99, normalized = TRUE, method = "parallel", nCores = 11)
# 
# indiv_saTBC_alpha0.01 <- get_sa_temporalBC(hypergraphList = incidMatList, timeStamps = timeStamps, smax = 25, windowLength = 2, 
#                                            focalTimeStamp = 2, alpha = 0.01, normalized = TRUE, method = "parallel", nCores = 11)

group_saTBC_alpha0.5 <- get_sa_temporalBC(hypergraphList = incidMatList, timeStamps = timeStamps, smax = 11, windowLength = 2, 
                                          focalTimeStamp = 2, alpha = 0.5, normalized = TRUE, method = "parallel", nCores = 20)

write.csv(group_saTBC_alpha0.5[[1]], file = file.path(sTBC_results, sprintf("sTBC_%s_alpha_0_dot_5.csv", run_ID)))

group_saTBC_alpha0.99 <- get_sa_temporalBC(hypergraphList = incidMatList, timeStamps = timeStamps, smax = 11, windowLength = 2, 
                                           focalTimeStamp = 2, alpha = 0.99, normalized = TRUE, method = "parallel", nCores = 20)

write.csv(group_saTBC_alpha0.99[[1]], file = file.path(sTBC_results, sprintf("sTBC_%s_alpha_0_dot_99.csv", run_ID)))

group_saTBC_alpha0.01 <- get_sa_temporalBC(hypergraphList = incidMatList, timeStamps = timeStamps, smax = 11, windowLength = 2, 
                                           focalTimeStamp = 2, alpha = 0.01, normalized = TRUE, method = "parallel", nCores = 20)

write.csv(group_saTBC_alpha0.01[[1]], file = file.path(sTBC_results, sprintf("sTBC_%s_alpha_0_dot_01.csv", run_ID)))

# write.csv(indiv_saTBC_alpha0.5[[1]], "C://Users/matth/Desktop/point5Dat2.csv")
# write.csv(indiv_saTBC_alpha0.99[[1]], "C://Users/matth/Desktop/point99Dat2.csv")
# write.csv(indiv_saTBC_alpha0.01[[1]], "C://Users/matth/Desktop/point01Dat2.csv")

# indiv_saTBC_alpha0.01 <- read.csv("C://Users/matth/Desktop/point01Dat2.csv", header = TRUE)
# indiv_saTBC_alpha0.5 <- read.csv("C://Users/matth/Desktop/point5Dat2.csv", header = TRUE)
# indiv_saTBC_alpha0.99 <- read.csv("C://Users/matth/Desktop/point99Dat2.csv", header = TRUE)
# 
# indivMat_T1 <- t(incidMatList[[1]])
# indivMat_T2 <- t(incidMatList[[2]])
# indiv_edgeSize_T1 <- sapply(1:ncol(indivMat_T1), function(x) sum(indivMat_T1[,x] > 0))
# indiv_edgeSize_T2 <- sapply(1:ncol(indivMat_T2), function(x) sum(indivMat_T2[,x] > 0))
# 
# edgeSizes <- sapply(1:length(indiv_edgeSize_T1), function(x) max(indiv_edgeSize_T1[x], indiv_edgeSize_T2[x]))
# 
# indiv_mean_saTBC_0.01 <- sapply(1:nrow(indiv_saTBC_alpha0.01), function(x) sum(indiv_saTBC_alpha0.01[x,2:26])/(ifelse(edgeSizes[x] > 25, 25, edgeSizes[x])))
# indiv_mean_saTBC_0.5 <- sapply(1:nrow(indiv_saTBC_alpha0.5), function(x) sum(indiv_saTBC_alpha0.5[x,2:26])/(ifelse(edgeSizes[x] > 25, 25, edgeSizes[x])))
# indiv_mean_saTBC_0.99 <- sapply(1:nrow(indiv_saTBC_alpha0.99), function(x) sum(indiv_saTBC_alpha0.99[x,2:26])/(ifelse(edgeSizes[x] > 25, 25, edgeSizes[x])))
# 
# indiv_meanSaTBC_Data <- data.table("ID" = indiv_saTBC_alpha0.01[,1], "alphaLow" = indiv_mean_saTBC_0.01,
#                                    "alphaMed" = indiv_mean_saTBC_0.5, "alphaHigh" = indiv_mean_saTBC_0.99)
# 
# indiv_meanSaTBC_Data$alphaLowRank <- rank(-indiv_meanSaTBC_Data$alphaLow)
# indiv_meanSaTBC_Data$alphaMedRank <- rank(-indiv_meanSaTBC_Data$alphaMed)
# indiv_meanSaTBC_Data$alphaHighRank <- rank(-indiv_meanSaTBC_Data$alphaHigh)
# 
# 
# 
# 
# jaccard(indiv_meanSaTBC_Data[alphaLowRank <= 40]$ID, indiv_meanSaTBC_Data[alphaMedRank <= 40]$ID)
# jaccard(indiv_meanSaTBC_Data[alphaLowRank <= 40]$ID, indiv_meanSaTBC_Data[alphaHighRank <= 40]$ID)
# jaccard(indiv_meanSaTBC_Data[alphaMedRank <= 40]$ID, indiv_meanSaTBC_Data[alphaHighRank <= 40]$ID)
# ids_0.01vs0.5 <- intersect(indiv_meanSaTBC_Data[alphaLowRank <= 40]$ID, indiv_meanSaTBC_Data[alphaMedRank <= 40]$ID)
# ids_0.01vs0.99 <- intersect(indiv_meanSaTBC_Data[alphaLowRank <= 40]$ID, indiv_meanSaTBC_Data[alphaHighRank <= 40]$ID)
# ids_0.5vs0.99 <- intersect(indiv_meanSaTBC_Data[alphaMedRank <= 40]$ID, indiv_meanSaTBC_Data[alphaHighRank <= 40]$ID)
# 
# cor.test(indiv_meanSaTBC_Data[ID %in% ids_0.01vs0.5]$alphaLow, indiv_meanSaTBC_Data[ID %in% ids_0.01vs0.5]$alphaMed, method = "kendall")
# cor.test(indiv_meanSaTBC_Data[ID %in% ids_0.01vs0.99]$alphaLow, indiv_meanSaTBC_Data[ID %in% ids_0.01vs0.99]$alphaHigh, method = "kendall")
# cor.test(indiv_meanSaTBC_Data[ID %in% ids_0.5vs0.99]$alphaMed, indiv_meanSaTBC_Data[ID %in% ids_0.5vs0.99]$alphaHigh, method = "kendall")

