source("./src/mcmcLinkReconcile/rcppMCMCDev/01-loadDataAndSources.R")
source('./src/mcmcLinkReconcile/rcppMCMCDev/00-Rfunctions.R')

# Compute pairs, get Ms and Us, and get the groupings

ms = c(.9, .7, .85, .4, .5, .7, .9, .9, .8)
us = sapply(demographics, computeU)
n = 4
priorLinkProb = .9
testData = demographics[1:n,]
testData = data.table(keyby = 1:n, testData[1:n])
setkey(testData, keyby)




# Comparison of proposal ratio functions


sexBlock = function(activeRecord, targetSet){
  all(activeRecord$gender == targetSet$gender)
}


#system.time(tmp <- linkConfigMHCompare(sample(1:10, n, replace = TRUE), checkBlocking = NULL, testData, ms, us, priorLinkProb, mcmc = 1000, burnin = 0, reportEvery = 1000))

#plot(apply(tmp, 1, function(x) length(unique(x))))

# Rcpp for log posterior

library(Rcpp)

# Compare logProp
labels = 1:3
naComparisons = compareAndValidateBlocking(c(1, 2, 2, 4), c(1, 2, 4, 4), 3, testData)
comparisons = naComparisons
comparisons[is.na(comparisons)] = 0

sourceCpp("./src/mcmcLinkReconcile/rcppMCMCDev/00-rcppFunctions.cpp")

microbenchmark(baseR = logConfigPost(1:3, naComparisons[,-c(1, 2)], ms, us, priorLinkProb), rcpp = logPosterior(1:3, naComparisons[,-c(1, 2)], ms, us, priorLinkProb))




# Compare blocking validation function

labelsCurrent = c(1, 2, 2, 4)
labelsPropose = c(2, 2, 2, 4)
changedLabelIndex = 1
dataMat = as.matrix(testData)



subsetAndValidateBlocking(labelsCurrent, labelsPropose, changedLabelIndex, dataMat)

microbenchmark(baseR = compareAndValidateBlocking(labelsCurrent, labelsPropose, changedLabelIndex, testData), rcpp = subsetAndValidateBlocking(labelsCurrent, labelsPropose, changedLabelIndex, dataMat))


# compare proposal ratio



baseRComps = compareAndValidateBlocking(labelsCurrent, labelsPropose, changedLabelIndex, testData)
rcppComps = subsetAndValidateBlocking(labelsCurrent, labelsPropose, changedLabelIndex, dataMat)

logProposalRatio(labelsCurrent, labelsPropose, comps, ms, us, priorLinkProb, changedLabelIndex)
logProp(labelsCurrent, labelsPropose, baseRComps, ms, us, priorLinkProb, changedLabelIndex)


microbenchmark(R = logProp(labelsCurrent, labelsPropose, baseRComps, ms, us, priorLinkProb, changedLabelIndex), cpp = logProposalRatio(labelsCurrent, labelsPropose, comps, ms, us, priorLinkProb, changedLabelIndex))




