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
              nMcmc = 1000,
              reportInterval = 10)

#Run MCMC

mcmcResults = partitionedMCMC(allBonds, compVectors, controlParameters = params)

# Profile the mcmc output, putting it into bond lists. 
mcmcProfile = lapply(mcmcResults, profileMcmc)

# find most probable configuration and its probability
mpc = findMpc(mcmcProfile)

# Profile the mpc

mpcProf = profileBonds(trueLabels, mpc$mpc, recs)
mpcProf






