source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

diffDataList <- importCSVs(path = "")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_diffCurveData="Sim-diffCurveData_higherOrderContagion_Weights_GoGvsHyp_paramSweep_AgeStrats_NetMets"
sim_T50="Sim-T50_higherOrderContagion_Weights_GoGvsHyp_paramSweep_AgeStrats_NetMets"
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
                                 "meanSeedAge" = 0, 
                                 "meanSeedDegree" = 0,
                                 "meanSeedBetween" = 0,
                                 "meanSeedStrength" = 0,
                                 "meanSeedsiD" = 0,
                                 "meanSeedsiBC" = 0,
                                 "meanSeedSED" = 0)

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
  diffusionCurveData$meanSeedDegree <- mean(sapply(initDemons, function(x) sum(diffDataTemp[which(diffDataTemp$ID == x),]$degree > diffDataTemp$degree)/99))
  diffusionCurveData$meanSeedBetween <- mean(sapply(initDemons, function(x) sum(diffDataTemp[which(diffDataTemp$ID == x),]$betweenness > diffDataTemp$betweenness)/99))
  diffusionCurveData$meanSeedStrength <- mean(sapply(initDemons, function(x) sum(diffDataTemp[which(diffDataTemp$ID == x),]$strength > diffDataTemp$strength)/99))
  diffusionCurveData$meanSeedsiD <- mean(sapply(initDemons, function(x) sum(diffDataTemp[which(diffDataTemp$ID == x),]$siD < diffDataTemp$siD)/99))
  diffusionCurveData$meanSeedsiBC <- mean(sapply(initDemons, function(x) sum(diffDataTemp[which(diffDataTemp$ID == x),]$siD < diffDataTemp$siBC)/99))
  diffusionCurveData$meanSeedSED <- mean(sapply(initDemons, function(x) sum(diffDataTemp[which(diffDataTemp$ID == x),]$subEdgeDens > diffDataTemp$subEdgeDens)/99))
  
  

  diffusionCurveData_T50 <- data.frame(
    "uniqueID" = diffusionCurveData$uniqueID[1],
    "topologyID" = diffusionCurveData$topologyID[1],
    "ageBias" = diffusionCurveData$ageBias[1],
    "selectGrad" = diffusionCurveData$selectGrad[1],
    "socialReinforcement" = diffusionCurveData$socialReinforcement[1],
    "seedStrategy" = diffusionCurveData$seedStrategy[1],
    "groupEffect" = diffusionCurveData$groupEffect[1],
    "T50" = min(diffusionCurveData[which(diffusionCurveData$propInformed >= 0.5),]$timeStep), 
    "meanSeedAge" = diffusionCurveData$meanSeedAge[1],
    "meanSeedDegree" = diffusionCurveData$meanSeedDegree[1],
    "meanSeedBetween" = diffusionCurveData$meanSeedBetween[1],
    "meanSeedStrength" = diffusionCurveData$meanSeedStrength[1],
    "meanSeedsiD" = diffusionCurveData$meanSeedsiD[1],
    "meanSeedsiBC" = diffusionCurveData$meanSeedsiBC[1],
    "meanSeedSED" = diffusionCurveData$meanSeedSED[1]
  )
  
  write.csv(diffusionCurveData, file = file.path(sim_diffCurveData, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(diffusionCurveData_T50, file = file.path(sim_T50, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  
}