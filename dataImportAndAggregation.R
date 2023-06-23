
importCSVs <- function() {
  # list all csv files from the current directory
  list.files(pattern=".csv$")
  # create a list from these files
  list.filenames<-list.files(pattern=".csv$")
  # create an empty list that will serve as a container to receive the incoming files
  list.data<-vector("list", length(list.filenames))
  
  # create a loop to read in your data
  for (i in 1:length(list.filenames)) {
    list.data[[i]]<-read.csv(list.filenames[i], header=TRUE)
    print(i)
    flush.console()
  }
  
  # add the names of your data to the list
  names(list.data)<-list.filenames
  
  return(list.data)
  
}

# aggregateData <- function(dataSets , indexAdjust) {
#   
#   summaryData <- data.frame(matrix(ncol = 1 + dim(dataSets[[1]])[2], nrow = dim(dataSets[[1]])[1] * length(dataSets)))
#   colnames(summaryData) <- c("uniqueID", colnames(list.data[[1]])[-1])
#   
#   for(i in 1:length(dataSets)) {
#     startIndex <- which(is.na(summaryData$uniqueID))[1]
#     endIndex <- dim(dataSets[[i]])[1] + startIndex - 1
#     summaryData[startIndex:endIndex,]$uniqueID <- i + indexAdjust
#     summaryData[startIndex:endIndex,2:dim(dataSets[[1]])[2]] <- dataSets[[i]][,2:dim(dataSets[[1]])[2]]
#   }
#   return(summaryData)
# }
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-details_simpleContagion_randomInherit/")
# list.data <- importCSVs()
# cSiR.Data <- aggregateData(dataSets = list.data, indexAdjust = 0)
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-details_simpleContagion_parentalInherit/")
# list.data <- importCSVs()
# cSiP.Data <- aggregateData(dataSets = list.data, indexAdjust = 250)
# 
# simpleContagionIndivData <- rbind(cSiR.Data, cSiP.Data)
# simpleContagionIndivData <- simpleContagionIndivData[,-43]
# 
# write.csv(simpleContagionIndivData, "C://Users/matth/Desktop/simpleContagionIndivData.csv")
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-details_numericContagion_randomInherit/")
# list.data <- importCSVs()
# cNCiR.Data <- aggregateData(dataSets = list.data, indexAdjust = 0)
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-details_numericContagion_parentalInherit/")
# list.data <- importCSVs()
# cNCiP.Data <- aggregateData(dataSets = list.data, indexAdjust = 250)
# 
# numericContagionIndivData <- rbind(cNCiR.Data, cNCiP.Data)
# numericContagionIndivData <- numericContagionIndivData[,-43]
# 
# write.csv(numericContagionIndivData, "C://Users/matth/Desktop/numericContagionIndivData.csv")
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-details_fractionalContagion_randomInherit/")
# list.data <- importCSVs()
# cFCiR.Data <- aggregateData(dataSets = list.data, indexAdjust = 0)
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-details_fractionalContagion_parentalInherit/")
# list.data <- importCSVs()
# cFCiP.Data <- aggregateData(dataSets = list.data, indexAdjust = 250)
# 
# fractionalContagionIndivData <- rbind(cFCiR.Data, cFCiP.Data)
# fractionalContagionIndivData <- fractionalContagionIndivData[,-43]
# 
# write.csv(fractionalContagionIndivData, "C://Users/matth/Desktop/fractionalContagionIndivData.csv")
# 
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-summaries_simpleContagion_parentalInherit/")
# cSiR.Data <- aggregateData(dataSets = list.data, indexAdjust = 0)
# cSiP.Data <- aggregateData(dataSets = list.data, indexAdjust = 250)
# 
# simpleContagionGlobalData <- rbind(cSiR.Data, cSiP.Data)
# simpleContagionGlobalData <- simpleContagionGlobalData[,-11]
# 
# write.csv(simpleContagionGlobalData, "C://Users/matth/Desktop/simpleContagionGlobalData.csv")
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-summaries_numericContagion_randomInherit/")
# list.data <- importCSVs()
# cNCiR.Data <- aggregateData(dataSets = list.data, indexAdjust = 0)
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-summaries_numericContagion_parentalInherit/")
# list.data <- importCSVs()
# cNCiP.Data <- aggregateData(dataSets = list.data, indexAdjust = 250)
# 
# numericContagionGlobalData <- rbind(cNCiR.Data, cNCiP.Data)
# numericContagionGlobalData <- numericContagionGlobalData[,-11]
# 
# write.csv(numericContagionGlobalData, "C://Users/matth/Desktop/numericContagionGlobalData.csv")
# 
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-summaries_fractionalContagion_randomInherit/")
# list.data <- importCSVs()
# cFCiR.Data <- aggregateData(dataSets = list.data, indexAdjust = 0)
# 
# setwd("C://Users/matth/Desktop/Simulation Results/Sim-summaries_fractionalContagion_parentalInherit/")
# list.data <- importCSVs()
# cFCiP.Data <- aggregateData(dataSets = list.data, indexAdjust = 250)
# 
# fractionalContagionGlobalData <- rbind(cFCiR.Data, cFCiP.Data)
# fractionalContagionGlobalData <- fractionalContagionGlobalData[,-11]
# 
# write.csv(fractionalContagionGlobalData, "C://Users/matth/Desktop/fractionalContagionGlobalData.csv")
