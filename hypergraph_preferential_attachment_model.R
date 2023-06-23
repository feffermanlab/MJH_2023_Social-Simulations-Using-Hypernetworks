source("Hypernetwork_functions.R")

generate_pa_hypergraph <- function(newVertexProb, initHypergraph, hyperEdgeSizes, targetSize, method, selectivity) {
  repeat{
    if(runif(1, min = 0, max = 1) <= newVertexProb) {
      newVertex <- max(unlist(initHypergraph)) + 1
      newEdgeSize <- sample(hyperEdgeSizes, size = 1, replace = TRUE)
      if(newEdgeSize > length(unique(unlist(initHypergraph)))) {
        newHyperEdge <- c(newVertex, unique(unlist(initHypergraph)))
        initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
      } else{
        degreeTemp <- sapply(sort(unique(unlist(initHypergraph))), function(x) sum(unlist(initHypergraph) == x))
        newHyperEdge <- c(newVertex, 
                          sample(sort(unique(unlist(initHypergraph))), size = (newEdgeSize - 1), replace = FALSE, prob = degreeTemp^selectivity))
        initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
      }
    } else {
      newEdgeSize <- sample(hyperEdgeSizes, size = 1, replace = TRUE)
      if(newEdgeSize > length(unique(unlist(initHypergraph)))) {
        newHyperEdge <- c(unique(unlist(initHypergraph)))
        initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
      } else {
        degreeTemp <- sapply(sort(unique(unlist(initHypergraph))), function(x) sum(unlist(initHypergraph) == x))
        newHyperEdge <- sample(sort(unique(unlist(initHypergraph))), size = newEdgeSize, 
                               replace = FALSE, prob = degreeTemp^selectivity)
        initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
      }
    }
    if(length(unique(unlist(initHypergraph))) >= max(hyperEdgeSizes)) {
      break
    }
  }
  # initVertices <- c(1:max(hyperEdgeSizes))
  # repeat{
  #   newEdgeSize <- sample(hyperEdgeSizes, size = 1, replace = TRUE)
  #   newHyperEdge <- sample(initVertices, size = newEdgeSize, replace = FALSE)
  #   initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
  #   if(length(unique(unlist(initHypergraph))) == length(initVertices)) {
  #     break
  #   }
  # }
  repeat{
    if(runif(1, min = 0, max = 1) <= newVertexProb) {
      newVertex <- max(unlist(initHypergraph)) + 1
      newEdgeSize <- sample(hyperEdgeSizes, size = 1, replace = TRUE)
      degreeTemp <- sapply(sort(unique(unlist(initHypergraph))), function(x) sum(unlist(initHypergraph) == x))
      newHyperEdge <- c(newVertex, 
                        sample(sort(unique(unlist(initHypergraph))), 
                               size = (newEdgeSize - 1), replace = FALSE, 
                               prob = degreeTemp^selectivity))
      initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
    } else{
      newEdgeSize <- sample(hyperEdgeSizes, size = 1, replace = TRUE)
      degreeTemp <- sapply(sort(unique(unlist(initHypergraph))), function(x) sum(unlist(initHypergraph) == x))
      newHyperEdge <- sample(sort(unique(unlist(initHypergraph))), size = newEdgeSize, 
                             replace = FALSE, prob = degreeTemp^selectivity)
      initHypergraph[[which(lengths(initHypergraph) == 0)[1]]] <- newHyperEdge
    }
    if(length(unique(unlist(initHypergraph))) >= targetSize) break
  }
  if(method == "incidence") {
    initHypergraph <- initHypergraph[-which(sapply(initHypergraph, is.null))]
    im <- get_incidence_matrix(hyperNetwork = initHypergraph, vertices = c(1:targetSize))
    return(im)
  } else{
    initHypergraph <- initHypergraph[-which(sapply(initHypergraph, is.null))]
  return(initHypergraph)
  }
}

################################
################################

hyperEdgeSizeDist <- read.csv("C://Users/matth/Desktop/hyperEdgeSizeDistribution.csv")
hyperEdgeSizes = hyperEdgeSizeDist[,2]

setwd("C://Users/matth/Desktop/Preferential Attachment Topologies/")

set.seed(5222023)
gL <- vector("list", 10000)
gL[[1]] <- 1
initHypergraph = gL

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-incidMat_diffusionTopologies"
if(!file.exists(sim_details)) dir.create(sim_details)

for(s in c(0,1,5)){
for(i in 1:100){
pa <- generate_pa_hypergraph(newVertexProb = 0.25, initHypergraph = gL, hyperEdgeSizes = hyperEdgeSizes, 
                              targetSize = 100, method = "incidence", selectivity = s)
write.csv(pa, file = file.path(sim_details, sprintf("incidMats_%s_%.01i_%.02i.csv", run_ID, s, i)))
}
}
