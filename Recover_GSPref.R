source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R")
library(foreach)
library(doParallel)
library(stringr)

registerDoParallel(cores = 20)


simData <- importCSVs(path = "~/scratch/SA_HyperNets/Run7/Sim-details_higherOrderContagion_Weights_GoGvsHyp_paramSweep_AgeStrats/")
popData <- importCSVs(path = "~/scratch/SA_HyperNets/Run7/Sim-livingPopData_diffusionTopologies_paramSweep/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details = "Sim-details_higherOrderContagion_Weights_GoGvsHyp_paramSweep_AgeStrats_withGSPref"
if(!file.exists(sim_details)) dir.create(sim_details)

foreach(s = 1:length(simData)) %dopar% {
  
  simDataTemp <- simData[[s]]
  
  simIDTemp <- simDataTemp$simID[1]
  
  popDataTemp <- popData[[simIDTemp]]
  
  if(sum(simDataTemp$ID == popDataTemp$ID) == length(simDataTemp$ID)) {
  simDataTemp$GSPref <- popDataTemp$GSPref
  } else{
    break
  }
  
  v = as.numeric(substring(names(simData[s]), 26, 28))
  g = as.numeric(substring(names(simData[s]), 30, 30))
  s_old = as.numeric(substring(names(simData[s]), 32, 32))
  
  write.csv(simDataTemp, file = file.path(sim_details, sprintf("simData_%s_%.1f_%.01i_%.01i_%.02i.csv", run_ID, v, g, s_old, simIDTemp)))
  
}