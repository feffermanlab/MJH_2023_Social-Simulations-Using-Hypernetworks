seedStrategy_highestDegree <- function(data, numSeeds) {
  data$degreeRank <- rank(-data$degree)
  newData <- data %>% arrange(data$degreeRank)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestBetweenness <- function(data, numSeeds) {
  data$betweennessRank <- rank(-data$betweenness)
  newData <- data %>% arrange(data$betweennessRank)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestharmCent <- function(data, numSeeds) {
  data$harmCentralityRank <- rank(-data$harmCentrality)
  newData <- data %>% arrange(data$harmCentralityRank)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestsiD <- function(data, numSeeds) {
  newData <- data %>% arrange(data$siD)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestsiBC <- function(data, numSeeds) {
  newData <- data %>% arrange(data$siBC)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestsiHC <- function(data, numSeeds) {
  newData <- data %>% arrange(data$siHC)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestGroupsiD <- function(data, numSeeds) {
  potentialSeeds <- as.integer(names(which(focalIncidMat[,data$topGroupsiD[1]] > 0)))
  if(length(potentialSeeds) > numSeeds) {
    seeds <- sample(potentialSeeds, size = numSeeds, replace = FALSE)
  } else{
    seeds <- potentialSeeds
  }
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestGroupsiBC <- function(data, numSeeds) {
  potentialSeeds <- as.integer(names(which(focalIncidMat[,data$topGroupsiBC[1]] > 0)))
  if(length(potentialSeeds) > numSeeds) {
    seeds <- sample(potentialSeeds, size = numSeeds, replace = FALSE)
  } else{
    seeds <- potentialSeeds
  }
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestGroupsiHC <- function(data, numSeeds) {
  potentialSeeds <- as.integer(names(which(focalIncidMat[,data$topGroupsiHC[1]] > 0)))
  if(length(potentialSeeds) > numSeeds) {
    seeds <- sample(potentialSeeds, size = numSeeds, replace = FALSE)
  } else{
    seeds <- potentialSeeds
  }
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}

seedStrategy_highestSED <- function(data, numSeeds) {
  newData <- data %>% arrange(-data$subEdgeDens)
  seeds <- newData$ID[1:numSeeds]
  data[which(data$ID %in% seeds),]$knowledgeState <- 1
  data[which(data$ID %in% seeds),]$initDemons <- 1
  return(data)
}
