# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
library(RecordLinkage)
library(R.utils)

sourceDirectory("./probLinkMetropolis/R/", modifiedOnly = FALSE)
sourceCpp("./probLinkMetropolis/src/00-rcppFunctions.cpp")

# Load Data 10000
data(RLdata10000)
recs = fread("./data/RLdata10000", colClasses = "character", na.strings = "")
recs = data.table(irow = 1:nrow(recs), recs[,c(1, 3, 5, 6, 7)])
trueLabels = identity.RLdata10000
# Prepare comparison vectors
allBonds = fread("./data/bondFiles/nineRuleBonds.csv")

# Load Data 500
# data(RLdata500)
# recs = fread("./data/RLdata500.csv", colClasses = "character", na.strings = "")
# recs = data.table(irow = 1:nrow(recs), recs[,c(1, 3, 5, 6, 7)])
# trueLabels = identity.RLdata500
# allBonds = fread("./data/bondFiles/rl500WithBlocks.csv")

compVectors = bondsToComparisonVectors(allBonds, recs)

# Set parameters

params = list(ms = c(.7, .9, .9, .9, .8),
              us = sapply(recs[,-1], computeU),
              priorLinkProb = 10^(-4),
              nMcmc = 100,
              reportInterval = 2000)

#Run MCMC

mcmcResults = partitionedMCMC(allBonds, compVectors, controlParameters = params)

# Profile the mcmc output, putting it into bond lists. 
mcmcProfile = lapply(mcmcResults, profileMcmc)

# find most probable configuration and its probability
mpc = findMpc(mcmcProfile)

# Profile the mpc

mpcProf = profileBonds(trueLabels, mpc$mpc, recs)
mpcProf

#
#
#

# Compute gamma from the full data set

trueM = computeM(recs, trueLabels)

# Run MCMC with Gibbs
gammaTotal = sapply(recs[,-1], function(col) sum(choose(table(col), 2)))
mcmcResultsGibbs <- partitionedMCMCGibbs(bondList = allBonds, comparisonVectors = compVectors, controlParameters = params, gammaTotal = gammaTotal, nc2 = choose(nrow(recs), 2))

# Trace plots
plot(mcmcResultsGibbs$alpha, type = 'l', ylab = "alpha")
matplot(mcmcResultsGibbs$m, ylab = "m", type = 'l', ylim = c(0, 1))
abline(h = trueM)
trueM
colMeans(mcmcResultsGibbs$m)
matplot(mcmcResultsGibbs$u, ylab = "u", type = 'l')

# Profile the mcmc output, putting it into bond lists. 
mcmcProfileGibbs = lapply(lapply(mcmcResultsGibbs$mcmc, function(x) x[-c(1:20),]), profileMcmc)

# find most probable configuration and its probability
mpcGibbs = findMpc(mcmcProfileGibbs)

# Profile the mpc

mpcProfGibbs = profileBonds(trueLabels, mpcGibbs$mpc, recs)
mpcProfGibbs

# Explore

comps = componentSizesFromBonds(allBonds)
noSingletons = lapply(which(comps$csize > 1), function(x) which(comps$membership == x))

rec = 2461
recComp = which(sapply(noSingletons, function(x) rec %in% x))
noSingletons[[recComp]]
mcmcProfileGibbs[[recComp]]
activeComps = compVectors[id1 %in% noSingletons[[recComp]]]
recs[irow %in% unlist(activeComps[,.(id1, id2)])]

