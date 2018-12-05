# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
library(RecordLinkage)

source("./probLinkMetropolis/R/functions.R")
sourceCpp("./probLinkMetropolis/src/00-rcppFunctions.cpp")
ari = mclust::adjustedRandIndex

# Load Data

data(RLdata500)
data500 = data.table(irow = 1:500, RLdata500[,c(1, 3, 5, 6, 7)]); rm(RLdata500)
trueLabels = identity.RLdata500

# Prepare comparisons

allBonds = disjunctionBlock(NULL, data500)

compVectors <- computeComparisonVectors(allBonds, as.matrix(data500))

# Set parameters

initialLabels = 1:nrow(data500)
ms = c(.7, .9, .9, .9, .8)
us = sapply(data500[,-1], computeU)
priorLinkProb = 1/uniqueN(data500[,-1])

t0 = Sys.time()
mhOut <- linkageMetropolis(initialLabels, compVectors, ms, us, priorLinkProb, mcmc = 10000, reportInterval = 100)
t1 = Sys.time() - t0
t1

system.time(aris <- apply(mhOut, 1, function(x) ari(x, trueLabels)))
plot(density(apply(mhOut, 1, uniqueN)))
plot(aris)

# Given a configuration and the truth, we want to compute some metrics to measure the performance, rand index is one

# This gives the number of bonds extracted by mcmc that were in the truth, and the number of true bonds that were captured by mcmc

profileBonds = function(trueLabels, mcmcLabels){
  
  nRuns = nrow(mcmcLabels)
  precision = numeric(nRuns)
  recall = numeric(nRuns)
  foundBondsList = vector("list", nRuns)
  falsePositiveList = vector("list", nRuns)
  
  trueBonds = data.table(t(apply(disjunctionBlock('trueLabels', data.table(irow = 1:length(trueLabels), trueLabels = trueLabels)), 1, sort)))
  setkey(trueBonds, V1, V2)
  
  for(i in 1:nRuns){
    
    mcmcBonds = data.table(t(apply(disjunctionBlock('mcmcLabels', data.table(irow = 1:length(trueLabels), mcmcLabels = mcmcLabels[i,])), 1, sort)))
    
    
    setkey(mcmcBonds, V1, V2)
    
    trueBonds[, foundByMcmc := FALSE][mcmcBonds, foundByMcmc := TRUE]
    mcmcBonds[, isTrueBond := FALSE][trueBonds, isTrueBond := TRUE]
    
    precision[i] = sum(mcmcBonds$isTrueBond)/nrow(mcmcBonds)
    recall[i] = sum(trueBonds$foundByMcmc)/nrow(trueBonds)
    foundBondsList[[i]] = trueBonds$foundByMcmc
    if(mean(mcmcBonds$isTrueBond) < 1) falsePositiveList[[i]] = data.table(mcmcBonds[isTrueBond == FALSE, list(V1, V2)], runIndex = i)
  }
  
  foundBonds = data.table(trueBonds[,list(V1, V2)], do.call('cbind', foundBondsList))
  colnames(foundBonds) = c("id1", "id2", paste0("Run", 1:nRuns))
  falsePositives = do.call('rbind', falsePositiveList)
  colnames(falsePositives) = c("id1", "id2", 'runIndex')
  
  out = list(recallMatrix = foundBonds, falsePositives = falsePositives, precision = precision, recall = recall)
  return(out)
}

rl500Profile = profileBonds(trueLabels, mhOut[1:1000,])
fpTable = rl500Profile$falsePositives[,.N, by = c("id1", "id2")][order(N)]
# 297 and 388 actually linked
cbind(trueLabels, data500)[c(297, 388, 410)]
table(rl500Profile$falsePositives[id2 == 410, 1:2])
rl500Profile$recallMatrix[id1 == 246 | id2 == 246]
