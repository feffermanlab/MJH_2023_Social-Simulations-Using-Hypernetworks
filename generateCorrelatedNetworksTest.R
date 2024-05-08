source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/Hypernetwork_functions.R", chdir = TRUE)
source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
library(gdata)
library(data.table)

generate_correlated_hypergraph <- function(hypergraph, rewireP) {
  numEdges <- ncol(hypergraph)
  rewire <- runif(numEdges, min = 0, max = 1)
  rewire <- ifelse(rewire <= rewireP, 1, 0)
  rewireList <- which(rewire == 1)
  if(length(rewireList) > 0) {
    rewireList <- resample(rewireList, length(rewireList))
    for(i in rewireList) {
      swapEdge <- sample(which(!(1:numEdges %in% i)),1)
      Hi <- resample(which(hypergraph[,i] > 0),1)
      Hx <- resample(which(hypergraph[,swapEdge] > 0),1)
      hypergraph[Hi,i] <- 0
      hypergraph[Hx,i] <- 1
    }
  }
  return(hypergraph)
}

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details_low="hyperGraph_Sequences_lowCorr_May24"
sim_details_high="hyperGraph_Sequences_highCorr_May24"
if(!file.exists(sim_details_low)) dir.create(sim_details_low)
if(!file.exists(sim_details_high)) dir.create(sim_details_high)

set.seed(05082024)

for(s in 1:100){

  H1 <- generate_ER_hypergraph(v = 500, m = 500, p = 0.01)

  for(t in 1:20) {

    fwrite(H1, file = file.path(sim_details_low, sprintf("simData_%s_%.01i_%.01i.csv", run_ID, s, t)))
    H1 <- generate_correlated_hypergraph(hypergraph = H1, rewireP = 0.5)

  }
}

set.seed(05082024)

for(s in 1:100){
  
  H1 <- generate_ER_hypergraph(v = 500, m = 500, p = 0.01)
  
  for(t in 1:20) {
    
    fwrite(H1, file = file.path(sim_details_high, sprintf("simData_%s_%.01i_%.01i.csv", run_ID, s, t)))
    H1 <- generate_correlated_hypergraph(hypergraph = H1, rewireP = 0.1)
    
  }
}