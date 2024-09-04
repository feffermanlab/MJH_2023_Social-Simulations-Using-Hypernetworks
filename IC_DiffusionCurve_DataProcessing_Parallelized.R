source("~/git/MJH_2023_Social-Simulations-Using-Hypernetworks/dataImportAndAggregation.R", chdir = TRUE)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)

#diffDataList <- importCSVs(path = "~/scratch/IC_Monograph/Sim-details_contagionResults_DirectMetrics_SFHH/")
#diffDataList <- importCSVs(path = "~/scratch/IC_Monograph/Sim-details_contagionResults_IndirectMetrics_SFHH/")
#diffDataList <- importCSVs(path = "~/scratch/IC_Monograph/Sim-details_contagionResults_DirectMetrics_InVS15/")
#diffDataList <- importCSVs(path = "~/scratch/IC_Monograph/Sim-details_contagionResults_IndirectMetrics_InVS15/")
diffDataList <- importCSVs(path = "~/scratch/IC_Monograph/Sim-details_contagionResults_DirectMetrics_LH10/")
#diffDataList <- importCSVs(path = "~/scratch/IC_Monograph/Sim-details_contagionResults_DirectMetrics_Thiers13/")

#Create folder in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
#sim_diffCurveData="Sim-diffCurveData_DirectMetrics_SFHH"
#sim_T50="Sim-T50_DirectMetrics_SFHH"
# sim_diffCurveData="Sim-diffCurveData_IndirectMetrics_SFHH"
# sim_T50="Sim-T50_IndirectMetrics_SFHH"
# sim_diffCurveData="Sim-diffCurveData_DirectMetrics_InVS15"
# sim_T50="Sim-T50_DirectMetrics_InVS15"
#sim_diffCurveData="Sim-diffCurveData_IndirectMetrics_InVS15"
#sim_T50="Sim-T50_IndirectMetrics_InVS15"
sim_diffCurveData="Sim-diffCurveData_DirectMetrics_LH10"
sim_T50="Sim-T50_DirectMetrics_LH10"
if(!file.exists(sim_diffCurveData)) dir.create(sim_diffCurveData)
if(!file.exists(sim_T50)) dir.create(sim_T50)

foreach(s = 1:length(diffDataList)) %dopar% {

diffDataTemp <- diffDataList[[s]]

maxTime <- max(diffDataTemp$acquisitionTime)

diffusionCurveData <- data.frame("uniqueID" = rep(s, maxTime),
                                 "topologyID" = 0,
                                 "socialReinforcement" = 0, 
                                 "seedStrategy" = 0,
                                 "timeStep" = 1:maxTime,
                                 "initDemons" = 0,
                                 "propInformed" = 0)

  diffusionCurveData$topologyID <- diffDataTemp$simID[1]
  diffusionCurveData$socialReinforcement <- diffDataTemp$socialReinforcement[1]
  diffusionCurveData$seedStrategy <- diffDataTemp$seedStrategy[1]
  initialNaive <- nrow(diffDataTemp[which(diffDataTemp$acquisitionTime > 0),])
  diffusionCurveData$initDemons <- nrow(diffDataTemp[which(diffDataTemp$initDemons == 1),])
  diffusionCurveData$propInformed <- sapply(seq(from = 1, to = maxTime), function(x) nrow(diffDataTemp[which(diffDataTemp$acquisitionTime <= x & diffDataTemp$acquisitionTime > 0),])/initialNaive)
  

  diffusionCurveData_T50 <- data.frame(
    "uniqueID" = diffusionCurveData$uniqueID[1],
    "topologyID" = diffusionCurveData$topologyID[1],
    "socialReinforcement" = diffusionCurveData$socialReinforcement[1],
    "seedStrategy" = diffusionCurveData$seedStrategy[1],
    "initDemons" = diffusionCurveData$initDemons[1],
    "T50" = min(diffusionCurveData[which(diffusionCurveData$propInformed >= 0.5),]$timeStep))
  
  write.csv(diffusionCurveData, file = file.path(sim_diffCurveData, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(diffusionCurveData_T50, file = file.path(sim_T50, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  
}