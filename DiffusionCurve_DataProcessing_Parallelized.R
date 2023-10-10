source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

#diffDataList <- importCSVs(path = "/home/mhasenja/scratch/SA_HyperNets/Run6/Sim-details_higherOrderContagion_Weights_GoGvsHyp_StratCompare/")
diffDataList <- importCSVs(path = "/home/mhasenja/scratch/SA_HyperNets/Run7/Sim-details_higherOrderContagion_Weights_GoGvsHyp_paramSweep/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_diffCurveData="Sim-diffCurveData_higherOrderContagion_Weights_GoGvsHyp_paramSweep"
sim_T50="Sim-T50_higherOrderContagion_Weights_GoGvsHyp_paramSweep"
if(!file.exists(sim_diffCurveData)) dir.create(sim_diffCurveData)
if(!file.exists(sim_T50)) dir.create(sim_T50)

foreach(s = 1:length(diffDataList)) %dopar% {

diffDataTemp <- diffDataList[[s]]

maxTime <- max(diffDataTemp$acquisitionTime)

diffusionCurveData <- data.frame("uniqueID" = rep(s, maxTime),
                                 "topologyID" = 0,
                                 "ageBias" = 0,
                                 "selectGrad" = 0,
                                 "socialReinforcement" = 0, 
                                 "seedStrategy" = 0,
                                 "groupEffect" = 0,
                                 "timeStep" = 1:maxTime,
                                 "propInformed" = 0,
                                 "meanSeedAge" = 0)

  diffusionCurveData$topologyID <- diffDataTemp$simID[1]
  diffusionCurveData$ageBias <- diffDataTemp$ageBias[1]
  diffusionCurveData$selectGrad <- diffDataTemp$selectGrad[1]
  diffusionCurveData$socialReinforcement <- diffDataTemp$socialReinforcement[1]
  diffusionCurveData$seedStrategy <- diffDataTemp$seedStrategy[1]
  diffusionCurveData$groupEffect <- diffDataTemp$groupEffect[1]
  initialNaive <- nrow(diffDataTemp[which(diffDataTemp$acquisitionTime > 0),])
  diffusionCurveData$propInformed <- sapply(seq(from = 1, to = maxTime), function(x) nrow(diffDataTemp[which(diffDataTemp$acquisitionTime <= x & diffDataTemp$acquisitionTime > 0),])/initialNaive)
  initDemons <- diffDataTemp[which(diffDataTemp$initDemons == 1),]$ID
  diffusionCurveData$meanSeedAge <- mean(sapply(initDemons, function(x) sum(x < diffDataTemp$ID)/99))

  diffusionCurveData_T50 <- data.frame(
    "uniqueID" = diffusionCurveData$uniqueID[1],
    "topologyID" = diffusionCurveData$topologyID[1],
    "ageBias" = diffusionCurveData$ageBias[1],
    "selectGrad" = diffusionCurveData$selectGrad[1],
    "socialReinforcement" = diffusionCurveData$socialReinforcement[1],
    "seedStrategy" = diffusionCurveData$seedStrategy[1],
    "groupEffect" = diffusionCurveData$groupEffect[1],
    "T50" = min(diffusionCurveData[which(diffusionCurveData$propInformed >= 0.5),]$timeStep), 
    "meanSeedAge" = diffusionCurveData$meanSeedAge[1]
  )
  
  write.csv(diffusionCurveData, file = file.path(sim_diffCurveData, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(diffusionCurveData_T50, file = file.path(sim_T50, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  
}