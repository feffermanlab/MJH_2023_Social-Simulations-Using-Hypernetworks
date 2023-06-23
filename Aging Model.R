#Could error be related to picking maxID as max of those alive?

#Load required functions
source("Hypernetwork_functions.R")
library(igraph)
library(stringr)

#Set working directory
setwd("C://Users/mjhas/OneDrive/Desktop/Social-Simulations-Using-Hypernetworks/")

#Set global parameters that do not change across simulations
#Parameters needed to generate initial hypergraph: are population size, number of hyperedges, and probability of hyperedge membership
N = 50
E = 25
p = 0.05

asocLearnProb = 0.05
maxT <- 200

numberSims = 10

#Can be "random" or "parental"
inheritance = "random"

#Create folders in which to store simulation results
#sim_details includes individual-level parameters at each time step for each simulation
#sim-summaries includes global-level statistics for each time step for each simulation run
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_details="Sim-details_numericContagion_randomInherit"
sim_summaries="Sim-summaries_numericContagion_randomInherit"
if(!file.exists(sim_details)) dir.create(sim_details)
if(!file.exists(sim_summaries)) dir.create(sim_summaries)

# ageBias = 0.2
# transmissionProbability = 0.05
# numericThreshold = 1
# fractionalThreshold = 0.33

#Set seed
set.seed(12012022)
set.seed(01252023)
set.seed(01282023)
set.seed(01302023)
set.seed(02012023)
set.seed(02062023)

ageBiasSet <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#transmissionStrengthSet <- c(0.05, 0.15, 0.25, 0.35, 0.45)
#transmissionStrengthSet <- c(0.1, 0.3, 0.5, 0.7, 0.9)
transmissionStrengthSet <- c(1, 3, 5, 7, 9)

#Set number of simulations to run
seedVect <- round(rnorm((numberSims * length(ageBiasSet) * length(transmissionStrengthSet)), mean = 500000, sd = 200000))

for(a in 1:length(ageBiasSet)) {
#for(a in 1){
  ageBias = ageBiasSet[a]
  for(m in 1:length(transmissionStrengthSet)) {
  #  for(m in 1:2) {
    transmissionProbability = transmissionStrengthSet[m]
    for(s in 1:numberSims) {
# ageBias = ageBiasSet[a]
# transmissionProbability = transmissionStrengthSet[m]
# ageBias = 0.1
# transmissionProbability = 0.05
# s = 1
      set.seed(seedVect[((a-1)*10)+m*s])
      
      initAges <- generate_age_structure(N, maxIter = 10000)
      popData <- create_population_data(N, initialAges = initAges, ageBias = ageBias)
      fullKnowledge <- generate_full_knowledge_set(numDomains = 5, numLevels = 5, complexity = 1)
      
      popNetwork <- generate_ER_hypergraph(v = N, m = E, p = p)
      popKnowledge <- generate_initial_pop_knowledge(popData = popData, knowledgeSet = fullKnowledge, asoc = asocLearnProb, simpleKnowledge = FALSE)
      
      #Next step is identify for each individual which individuals it has previously associated with
      #To do this, I need to transform the hypernetwork to its dual and extract its linegraph with s = 1
      #Then, an individual will be directly connected to its associates in the line graph
      dualPopNetwork <- get_dual_hypergraph(hypNet = popNetwork, popData = popData)
      
      #Produces s-line graph where the vertices correspond to the individuals
      dualLineGraph <- get_s_line_graph(hypNet = dualPopNetwork, size = 1)
      
      #Could try preallocating size of individual matrices in the array
      preferenceMatrices <- initialize_preference_matrix(N = N, maxT = maxT)
      
      individual_summary_data <- data.frame("simID" = s, 
                                            "TransModel" = "simple",
                                            "Ptrans" = transmissionProbability,
                                            "ageBias" = ageBias,
                                            "Inheritance" = inheritance,
                                            "ID" = 0, 
                                            "timeStep" = 0, 
                                            "Age" = 0,
                                            "Parent" = 0,
                                            "GSPref" = 0,
                                            "AgeBias" = 0,
                                            "NVariants" = 0,
                                            "NSimpleVariants" = 0,
                                            "N2Variants" = 0,
                                            "variantMax" = rep(length(fullKnowledge), N * maxT),
                                            "simpleVariantMax" = rep(sum(sapply(seq(1:length(fullKnowledge)), function(x) sum(sapply(fullKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == 1)))), 
                                                                     N * maxT),
                                            "C2VariantMax" = rep(sum(sapply(seq(1:length(fullKnowledge)), function(x) sum(sapply(fullKnowledge[[x]], function(a) str_count(a, substring(a,1,1)) == 2)))), 
                                                                 N * maxT),
                                            "s1BC" = 0,
                                            "s2BC" = 0,
                                            "s3BC" = 0,
                                            "s4BC" = 0,
                                            "s5BC" = 0,
                                            "s6BC" = 0,
                                            "s7BC" = 0,
                                            "s8BC" = 0,
                                            "s9BC" = 0,
                                            "s10BC" = 0,
                                            "siBC" = 0,
                                            "strength" = 0, 
                                            "wBC" = 0, 
                                            "s1D" = 0,
                                            "s2D" = 0,
                                            "s3D" = 0,
                                            "s4D" = 0,
                                            "s5D" = 0,
                                            "s6D" = 0,
                                            "s7D" = 0,
                                            "s8D" = 0,
                                            "s9D" = 0,
                                            "s10D" = 0,
                                            "siD" = 0)
      
      global_summary_data <- data.frame("simID" = s, 
                                        "TransModel" = "simple",
                                        "Ptrans" = transmissionProbability,
                                        "ageBias" = ageBias,
                                        "Inheritance" = inheritance,
                                        "timeStep" = seq(1:maxT), 
                                        "knowledgeRichness" = 0, 
                                        "knowledgeDiversity" = 0, 
                                        "knowledgeEvenness" = 0)
      for(t in 1:maxT){
        timeStep <- t
        
        currentPartners <- vector("list", 12 * nrow(popData[which(popData$Alive == "Y"),]))
        
        for(k in 1:12){
          partners <- select_partners(prefMatrix = preferenceMatrices, popData = popData, t = timeStep)
          newKnowledge <- complex_contagion_numeric_threshold(hyperEdges = partners, popKnowl = popKnowledge, 
                                                              popData = popData, threshold = transmissionProbability)
          # newKnowledge <- complex_contagion_fractional_threshold(hyperEdges = partners, popKnowl = popKnowledge, 
          #                                                        popData = popData, threshold = transmissionProbability)
            #simple_contagion(hyperEdges = partners, popKnowl = popKnowledge, popData = popData, probTrans = transmissionProbability)
          if(sum(!(sapply(newKnowledge, is.null))) > 0) {
            popKnowledge <- update_knowledge(popKnowl = popKnowledge, newKnowl = newKnowledge)
          }
          indexStart <- (length(currentPartners) - length(currentPartners[which(sapply(currentPartners, is.null))])) + 1
          indexStop <- indexStart + length(partners) - 1
          currentPartners[indexStart:indexStop] <- partners
          }
        currentPartners <- currentPartners[-which(sapply(currentPartners, is.null))]
        
        individual_summary_data[(((t - 1) * N) + 1):(N * t), 
                                6:41] <- record_individual_summary_data(popData = popData, popKnowledge = popKnowledge, knowledgeSet = fullKnowledge)
        global_summary_data[t,6:9] <- record_global_summary_data(popData = popData, popKnowledge = popKnowledge, timeStep = t)
        
        demoChanges <- identify_demographic_changes(popData = popData)
        
        preferenceMatrices[[timeStep + 1]] <- update_preference_matrix(prefMatrices = preferenceMatrices, popData = popData, timeStep = timeStep, 
                                                                       demoChanges = demoChanges, currentPartners = currentPartners, inheritance = inheritance)
        
        #Could creating popData anew here be slowing the loop down?
        popData <- update_population_data(popData = popData, demoChanges = demoChanges, ageBias = ageBias)
        if(length(demoChanges[[1]]) > 0) {
          newIndivsKnowl <- generate_knowledge_for_new_indivs(popData = popData, knowledgeSet = fullKnowledge, asoc = asocLearnProb, simpleKnowledge = TRUE)
          for(i in 1:length(which(popData$Age == 1))) {
            popKnowledge[[popData$ID[which(popData$Age == 1)][i]]] <- newIndivsKnowl[[i]]
          }
        }
        
        popKnowledge <- asocial_learning(popKnowledge = popKnowledge, popData = popData, knowledgeSet = fullKnowledge, simpleKnowledge = TRUE, asoc = asocLearnProb)
        
      }
      
      popData <- popData[which(popData$Age > 0),]
      popKnowledge <- popKnowledge[-which(sapply(popKnowledge, is.null))]
      
      write.csv(individual_summary_data, file = file.path(sim_details, sprintf("ILData_%s_%.2f_%.2f_%02i.csv", run_ID, ageBias, transmissionProbability, s)))
      write.csv(global_summary_data, file = file.path(sim_summaries, sprintf("GLData_%s_%.2f_%.2f_%02i.csv", run_ID, ageBias, transmissionProbability, s)))
      
    }
  }
}


# newKnowledgeCN <- complex_contagion_numeric_threshold(hyperEdges = currentPartners, popKnowl = popKnowledge, popData = popData, threshold = numericThreshold)
# newKnowledgeCF <- complex_contagion_fractional_threshold(hyperEdges = currentPartners, popKnowl = popKnowledge, popData = popData, threshold = fractionalThreshold)

# hyperEdgeKnowledge <- get_hyperedge_knowledge(hyperEdges = currentPartners, popKnowl = popKnowledge)
# 
# newKnowledge <- learn_new_knowledge(hyperEdgeKnowledge = hyperEdgeKnowledge, popData = popData, popKnowl = popKnowledge, hyperEdges = currentPartners, 
#                                     selectMethod = "total", kThreshold = 0.5)