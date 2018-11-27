source("./src/mcmcLinkReconcile/rcppMCMCDev/01-loadDataAndSources.R")
source('./src/mcmcLinkReconcile/rcppMCMCDev/00-Rfunctions.R')
library(Rcpp)
sourceCpp("./src/mcmcLinkReconcile/rcppMCMCDev/00-rcppFunctions.cpp")

# Compute pairs, get Ms and Us, and get the groupings

ms = c(.9, .7, .85, .4, .5, .7, .9, .9, .8)
us = sapply(demographics, computeU)
n = 4
priorLinkProb = .9
testData = demographics[1:n,]
testData = data.table(keyby = 1:n, testData[1:n])
setkey(testData, keyby)
dataMat = as.matrix(testData)



# Compute all valid comparison pairs according to the blocking criteria

sexBlock = function(x1, x2){
  index = which(names(x1) == 'gender')
  x1[index] == x2[index]
}

extractValidComparisons = function(data, blockingFn){
  comparisons = list()
  pairs = t(combn(1:nrow(data), 2))
  nc2 = nrow(pairs)
  pairsMatched = apply(pairs, 1, function(x){
    blockingFn(data[x[1],], data[x[2],])
  })
  return(pairs[pairsMatched,])
}

blockedComps = extractValidComparisons(dataMat, sexBlock)

# Pass this matrix along with the data to cpp, compute the comparison vectors

comps = computeComparisonVectors(blockedComps, dataMat)

# Now pass some labels to logProposalRatio. It should return 'reject' if the needed comparisons can't be found


labelsCurrent = c(2, 2, 2, 3)
labelsPropose = c(2, 2, 3, 3)
changedLabelIndex = 3
baseRComps = compareAndValidateBlocking(labelsCurrent, labelsPropose, changedLabelIndex, testData)


logProposalRatio(labelsCurrent, labelsPropose, comps, ms, us, priorLinkProb, changedLabelIndex)
logProp(labelsCurrent, labelsPropose, baseRComps, ms, us, priorLinkProb, changedLabelIndex)


#microbenchmark(R = logProp(labelsCurrent, labelsPropose, baseRComps, ms, us, priorLinkProb, changedLabelIndex), cpp = logProposalRatio(labelsCurrent, labelsPropose, comps, ms, us, priorLinkProb, changedLabelIndex))

## Test MH

#microbenchmark(r = linkConfigMH(1:4, baseRComps, ms, us, priorLinkProb, mcmc = 5000, reportEvery = 10000), cpp = linkConfigMHCpp(1:4, comps, ms, us, priorLinkProb, mcmc = 5000, reportEvery = 10000), times = 10)

sourceCpp("./src/mcmcLinkReconcile/rcppMCMCDev/00-rcppFunctions.cpp")

linkageMetropolis(c(1, 2, 2, 3), comps, ms, us, priorLinkProb, mcmc = 10, reportInterval = 10000)

tm = microbenchmark(fullrcpp = linkageMetropolis(c(1, 2, 2, 3), comps, ms, us, priorLinkProb, mcmc = 5000, reportInterval = 10000), partcpp = linkConfigMHCpp(1:4, comps, ms, us, priorLinkProb, mcmc = 5000, reportEvery = 10000), justR = linkConfigMH(1:4, baseRComps, ms, us, priorLinkProb, mcmc = 5000, reportEvery = 10000), times = 10)

autoplot.microbenchmark(tm)

