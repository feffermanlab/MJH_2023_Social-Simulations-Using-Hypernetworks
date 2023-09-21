
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
library(igraph)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_networkData="Sim-networkData_diffusionTopologies_paramSweep"
if(!file.exists(sim_networkData)) dir.create(sim_networkData)

popDataList <- importCSVs(path = "/home/mhasenja/scratch/SA_HyperNets/Run5/Sim-livingPopData_diffusionTopologies_paramSweep/")
incidMatList <- importCSVs(path = "/home/mhasenja/scratch/SA_HyperNets/Run5/Sim-incidMat_diffusionTopologies_paramSweep/")

popData <- data.frame("simID" = rep(0, 100),
                      "ID" = 0,
                      "ID.SimID" = 0,
                      "timeStep" = 0,
                      "ageBias" = 0,
                      "selectGrad" = 0,
                      "GSPref" = 0, 
                      "Parent" = 0,
                      "Age" = 0,
                      "BIAge" = 0,
                      "WIAge" = 0,
                      "Degree" = 0,
                      "Between" = 0,
                      "eVC" = 0,
                      "siDegree" = 0,
                      "siBetween" = 0,
                      "subEdgeDens" = 0)

foreach(i = 1:length(popDataList)) %dopar% {
  # startIndex = (i * 100) - 99
  # endIndex  = (i * 100)
  simIDTemp = ceiling(i/10)
  nameLength = nchar(names(popDataList)[[i]])
  timeTemp = as.integer(substring(names(popDataList)[[i]], nameLength - 6, nameLength - 4))
  #popData[startIndex:endIndex, "ID"] <- popDataList[[i]]$ID
  popData$ID <- popDataList[[i]]$ID
  popData$simID <- simIDTemp
  popData$ID.SimID <- paste(popData$ID, popData$simID, sep = ".")
  popData$timeStep <- timeTemp
  popData$ageBias <- popDataList[[i]]$AgeBias
  popData$selectGrad <- popDataList[[i]]$selectGradient
  popData$GSPref <- popDataList[[i]]$GSPref
  popData$Parent <- popDataList[[i]]$Parent
  popData$Age <- popDataList[[i]]$Age
  
  currentIDs <- incidMatList[[i]][,1]
  focalIncidMat <- as.matrix(incidMatList[[i]][,-1])
  focalIncidMat <- as.matrix(focalIncidMat[,-dim(focalIncidMat)[2]])
  row.names(focalIncidMat) <- currentIDs
  colnames(focalIncidMat) <- seq(from = 1, to = dim(focalIncidMat)[2])
  dualIncidMat <- t(focalIncidMat)
  GoGmatrix <- get_adjacency_matrix_from_incidence_matrix(I = focalIncidMat, 
                                                          V = currentIDs, 
                                                          weighted = TRUE)
  GoGgraph <- graph_from_adjacency_matrix(GoGmatrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  popData$Degree <- as.vector(degree(GoGgraph))
  popData$Between <- as.vector(betweenness(GoGgraph, 
                                                                   directed = FALSE, normalized = TRUE, 
                                                                   weights = 1/edge.attributes(GoGgraph)$weight))
  popData$eVC <- as.vector(eigen_centrality(GoGgraph, directed = FALSE, scale = TRUE)$vector)
  
  popData$siDegree <- get_s_degree(hypergraph = dualIncidMat, smax = 6, vertexNames = currentIDs, mode = "incidence")[[2]]
  popData$siBetween <- get_s_betweenness(hypergraph = dualIncidMat, smax = 6, vertexNames = currentIDs, mode = "incidence")[[2]]
  
  popData$subEdgeDens <- sapply(popData$ID, function(x) get_local_subedge_density(hypergraph = focalIncidMat, vertex = x))
  
  write.csv(popData, file = file.path(sim_networkData, sprintf("popData_%s_%.2f_%02i.csv", run_ID, simIDTemp, i)))
}
