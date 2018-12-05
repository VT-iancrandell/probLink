# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
source("./probLinkMetropolis/R/functions.R")
sourceCpp("./probLinkMetropolis/src/00-rcppFunctions.cpp")
demographics = fread("./data/cleanedDHSData.csv", na.strings = '')
setkey(demographics, irow)

blockingVars = c("last_name")

tmp = demographics[1:1000]
allBonds = disjunctionBlock(blockingVars, tmp)

compVectors <- computeComparisonVectors(allBonds, as.matrix(tmp))



