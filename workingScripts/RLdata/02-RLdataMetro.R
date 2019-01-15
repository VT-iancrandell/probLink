# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
library(RecordLinkage)
library(R.utils)

sourceDirectory("./probLinkMetropolis/R/", modifiedOnly = FALSE)
sourceCpp("./probLinkMetropolis/src/00-rcppFunctions.cpp")

# Load Data
data(RLdata500)
recs = fread("./data/RLdata500.csv", colClasses = "character", na.strings = "")
recs = data.table(irow = 1:500, recs[,c(1, 3, 5, 6, 7)])
trueLabels = identity.RLdata500

# Prepare comparison vectors

allBonds = fread("./data/bondFiles/rl500WithBlocks.csv")
compVectors = bondsToComparisonVectors(allBonds, recs)

# Set parameters

params = list(ms = c(.7, .9, .9, .9, .8),
              us = sapply(recs[,-1], computeU),
              priorLinkProb = 10^(-4),
              nMcmc = 100,
              reportInterval = 10)

#Run MCMC

mcmcResults = partitionedMCMC(allBonds, compVectors, controlParameters = params)

# find most probable configuration and its probability

mcmcComponent = mcmcResults[[29]]

# This is written to run on a single component at a time
profileMcmc = function(componentMcmc){
  
  mcmcBonds = lapply(1:nrow(componentMcmc), function(x){
    row = sort(componentMcmc[x,])
    out = do.call(rbind, componentLabelsToBondsCharacter(row, as.numeric(names(row)), uniqueN(row)))
    out = out[order(out[,1]),]
  })
  
  configurationIndices = numeric(length(mcmcBonds))
  configurationIndices[1] = 1
  
  configurations = vector("list")
  configurationIndex = 1
  configurations[[configurationIndex]] = mcmcBonds[[1]]
  
  for(i in 2:length(mcmcBonds)){
    
    addNewPrototype = TRUE
    
    for(j in 1:configurationIndex){
      if(identical(mcmcBonds[[i]], configurations[[j]])){
        addNewPrototype = FALSE
        break
      }
    }
    if(addNewPrototype){
      configurationIndex = configurationIndex + 1
      configurations[[configurationIndex]] = mcmcBonds[[i]]
      configurationIndices[i] = configurationIndex
    }else{
      configurationIndices[i] = j
    }
  }
  
  configurationFrequencies = data.table(configurationIndex = 1:configurationIndex, frequency = tabulate(configurationIndices)/nrow(componentMcmc))
  return(list(configurations = configurations, whichPrototype = configurationIndices, configurationFrequencies = configurationFrequencies))
}

prof = profileMcmc(mcmcComponent)
prof




head(tmp)
head(mcmcComponent)


mpc = findMPC(mcmcResults)

# Profile the mpc

mpcProf = profileBonds(trueLabels, mpc$componentID, recs)
mpcProf






