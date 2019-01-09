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
recs = data.table(irow = 1:nrow(recs), recs[,c(1, 3, 5, 6, 7)])
trueLabels = identity.RLdata500
n = nrow(recs)

# Prepare comparison vectors

allBonds = fread("./data/bondFiles/rl500NoBlocks.csv")

compVectors = bondsToComparisonVectors(allBonds, recs)

# Set parameters

params = list(ms = c(.7, .9, .9, .9, .8),
              us = sapply(recs[,-1], computeU),
              priorLinkProb = 1/500,
              nMcmc = 1000)


#Run MCMC
# State = Z vector and M vector, iterate through one than the other.
# This is for one component at a time
# Need to keep track of links as the algorithm progresses, so update link tables at every step
# Each label (and set of obs with that label) have a 

oldLabels = 1:nrow(recs)



# All together

mcmcResultsGibbs = linkageMetropolisWithGibbs(oldLabels, as.matrix(compVectors), params$ms, params$us, params$priorLinkProb, mcmc = 1000, reportInterval = 100)
names(mcmcResultsGibbs) = c("labels", "alpha", "m", "u")
colnames(mcmcResultsGibbs$labels) = 1:500
matplot(mcmcResultsGibbs$alpha, type = 'l', ylim = c(0, 1))
matplot(mcmcResultsGibbs$m, type = 'l', ylim = c(0, 1))
matplot(mcmcResultsGibbs$u, type = 'l', ylim = c(0, 1))

profileGibbs = profileMcmcBonds(trueLabels, mcmcResultsGibbs$labels)
#mpcGibbs = findMPC(list(mcmcResultsGibbs$labels))

rowMeans(profileGibbs$recallMatrix[,-(1:2)])

plot(profileGibbs$precision, col = 'red', type = 'l', ylim = c(0 ,1))
lines(profileGibbs$recall, col = 'blue', type = 'l')
legend("bottomright", legend = c("Precision", "Recall"), fill = c("red", "blue"), bty = 'n')
mpc = findMPC(list(mcmcResultsGibbs$labels[999:1000,]))

mpcProf = profileBonds(trueLabels, mpc$componentID, recs)








