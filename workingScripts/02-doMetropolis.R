# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
source("./workingScripts/01-makeComparisonVectors.R")
demographics = fread("./data/cleanedDHSData.csv", na.strings = '')
setkey(demographics, irow)


initialLabels = 1:nrow(tmp)
ms = c(.9, .7, .85, .4, .5, .7, .9, .9, .8)
us = sapply(tmp[,-1], computeU)
priorLinkProb = 1/uniqueN(tmp)

system.time(mhOut <- linkageMetropolis(initialLabels, compVectors, ms, us, priorLinkProb, mcmc = 2000, reportInterval = 100))
