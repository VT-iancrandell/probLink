# This script incorporates the blocking rules and possible deduped database assumption to create the comparison vectors for linkage
library(Rcpp)
library(igraph)
library(data.table)
library(RecordLinkage)
library(R.utils)

sourceDirectory("./probLinkMetropolis/R/", modifiedOnly = FALSE)
sourceCpp("./probLinkMetropolis/src/00-rcppFunctions.cpp")

# Load Data
data(RLdata10000)
recs = fread("./data/RLdata10000", colClasses = "character", na.strings = "")
recs = data.table(irow = 1:10000, recs[,c(1, 3, 5, 6, 7)])
trueLabels = identity.RLdata10000

# Prepare comparison vectors

allBonds = fread("./data/bondFiles/nineRuleBonds.csv")

compVectors = bondsToComparisonVectors(allBonds, recs)

# Set parameters

params = list(ms = c(.7, .9, .9, .9, .8),
              us = sapply(recs[,-1], computeU),
              priorLinkProb = 10^(-4),
              nMcmc = 100)


#Run MCMC

mcmcResults = partitionedMCMC(allBonds, compVectors, controlParameters = params)

# This gives the number of bonds extracted by mcmc that were in the truth, and the number of true bonds that were captured by mcmc
profile = profileMcmcBonds(trueLabels, mcmcResults)


plot(profile$precision, col = 'red', type = 'l', ylim = c(0 ,1))
lines(profile$recall, col = 'blue', type = 'l')
legend("bottomright", legend = c("Precision", "Recall"), fill = c("red", "blue"), bty = 'n')


# find most probable configuration and its probability

mpc = findMPC(mcmcResults)

# Profile the mpc steam

mpcProf = profileBonds(trueLabels, mpc$componentID, recs)
mpcProf

for(i in 1:5) print(recs[irow %in% mpcProf$falsePositives[i]])






