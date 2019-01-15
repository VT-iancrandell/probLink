# Wrapper for the MCMC function to execust on a list of partitioned elements
#
# WORKING
#
# 
# bondList = allBonds
# comparisonVectors = compVectors
# controlParameters = list(ms = c(.7, .9, .9, .9, .8),
#                          us = sapply(recs[,-1], computeU),
#                          priorLinkProb = 10^(-4),
#                          nMcmc = 100,
#                          reportInterval = 100)

partitionedMCMCGibbs = function(bondList, comparisonVectors, controlParameters){
  
  # Unpack arg list
  
  ms = controlParameters$ms
  us = controlParameters$us
  priorLinkProb = controlParameters$priorLinkProb
  nMcmc = controlParameters$nMcmc
  reportInterval = controlParameters$reportInterval
  
  comps = componentSizesFromBonds(bondList)
  noSingletons = lapply(which(comps$csize > 1), function(x) which(comps$membership == x))
  nComps = length(noSingletons)
  gammaTotal = colSums(comparisonVectors[,-c(1, 2)])
  
  
  mcmcOut = matrix(NA, nMcmc, sum(sapply(noSingletons, length)))
  alaphOut = numeric(nMcmc)
  mOut = matrix(NA, nMcmc, length(ms))
  uOut = matrix(NA, nMcmc, length(us))
  
  
  for(i in 1:nMcmc){
    
    partitionedComparisonVectors = lapply(1:nComps, function(x) {
      
      initialLabels = 1:length(noSingletons[[x]])
      activeComp = comparisonVectors[id1 %in% noSingletons[[x]]]
      relabeling = cbind(match(activeComp$id1, noSingletons[[x]]), match(activeComp$id2, noSingletons[[x]]))
      out = linkageMetropolis(initialLabels, as.matrix(cbind(relabeling, activeComp[,-c(1, 2)])), ms, us, priorLinkProb, mcmc = 1, reportInterval = reportInterval)
      colnames(out) = noSingletons[[x]]
      apply(out, 2, function(y) paste0(y, ".", x))
    })
    
    mcmcOut[i,] = unlist(partitionedComparisonVectors)
    
    bondList = do.call(rbind, 
                       lapply(1:nComps, 
                              function(x){
                                makeBonds(data.table(labels = partitionedComparisonVectors[[x]], 
                                                                irow = names(partitionedComparisonVectors[[x]])), "labels")
                              }
                       )
    )
    bondList = apply(bondList, 1, function(x) as.numeric(x))
    
    gibbsOut = gibbsStep(bondList, as.matrix(comparisonVectors), gammaTotal)
    alaphOut[i] = gibbsOut[[1]]
    mOut[i,] = gibbsOut[[2]]
    uOut[i,] = gibbsOut[[3]]
    
  }
  return(list(alpha = alphaOut, m = mOut, u = uOut, mcmc = mcmcOut))
}

