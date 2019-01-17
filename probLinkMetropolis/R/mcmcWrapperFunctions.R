# Run the no-gibbs version of the algorithm on a list

partitionedMCMC = function(bonds, comparisonVectors, controlParameters){
  
  # Unpack arg list
  
  ms = controlParameters$ms
  us = controlParameters$us
  priorLinkProb = controlParameters$priorLinkProb
  nMcmc = controlParameters$nMcmc
  reportInterval = controlParameters$reportInterval
  
  comps = componentSizesFromBonds(bonds)
  noSingletons = lapply(which(comps$csize > 1), function(x) which(comps$membership == x))
  nComps = length(noSingletons)
  
  mcmcLabels = lapply(1:nComps, function(x) {
    
    initialLabels = 1:length(noSingletons[[x]])
    activeComp = comparisonVectors[id1 %in% noSingletons[[x]]]
    relabeling = cbind(match(activeComp$id1, noSingletons[[x]]), match(activeComp$id2, noSingletons[[x]]))
    out = linkageMetropolis(initialLabels, as.matrix(cbind(relabeling, activeComp[,-c(1, 2)])), ms, us, priorLinkProb, mcmc = nMcmc, reportInterval = reportInterval)
    colnames(out) = noSingletons[[x]]
    apply(out, 2, function(y) paste0(y, ".", x))
  })
}



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
# gammaTotal = gammaTotal
# nc2 = choose(nrow(recs), 2)
partitionedMCMCGibbs = function(bondList, comparisonVectors, controlParameters, gammaTotal, nc2){
  
  # Unpack arg list
  
  ms = controlParameters$ms
  us = controlParameters$us
  priorLinkProb = controlParameters$priorLinkProb
  nMcmc = controlParameters$nMcmc
  reportInterval = controlParameters$reportInterval
  
  comps = componentSizesFromBonds(bondList)
  noSingletons = lapply(which(comps$csize > 1), function(x) which(comps$membership == x))
  nComps = length(noSingletons)
  if(missing(gammaTotal)) gammaTotal = colSums(comparisonVectors[,-c(1, 2)])
  
  
  mcmcOut = vector("list", nMcmc)
  alphaOut = numeric(nMcmc)
  mOut = matrix(NA, nMcmc, length(ms))
  uOut = matrix(NA, nMcmc, length(us))
  oldLabels = vector("list", nMcmc)
  
  for(i in 1:nMcmc){
    
    mcmcLabels = lapply(1:nComps, function(x) {
      
      if(i == 1){
        initialLabels = 1:length(noSingletons[[x]]) 
      }else{
        initialLabels = oldLabels[[x]]
      }
      activeComp = comparisonVectors[id1 %in% noSingletons[[x]]]
      relabeling = cbind(match(activeComp$id1, noSingletons[[x]]), match(activeComp$id2, noSingletons[[x]]))
      oldLabels[[x]] <<- linkageMetropolis(initialLabels, as.matrix(cbind(relabeling, activeComp[,-c(1, 2)])), ms, us, priorLinkProb, mcmc = 1, reportInterval = reportInterval)
      colnames(oldLabels[[x]]) = noSingletons[[x]]
      apply(oldLabels[[x]], 2, function(y) paste0(y, ".", x))
    })
    mcmcOut[[i]] = mcmcLabels
    mcmcBonds <- do.call(rbind,
                        lapply(1:nComps, function(x){
                          row = sort(mcmcLabels[[x]])
                          out = do.call(rbind, componentLabelsToBondsCharacter(row, as.numeric(names(row)), uniqueN(row)))
                          out = out[order(out[,1]),, drop = FALSE]
                        }))
    
    gibbsOut = gibbsStep(mcmcBonds, as.matrix(comparisonVectors), gammaTotal, nc2)
    priorLinkProb = alphaOut[i] = gibbsOut[[1]]
    ms = mOut[i,] = gibbsOut[[2]]
    us = uOut[i,] = gibbsOut[[3]]
    
  }
  
  mcmcOut = lapply(1:nComps, function(x){
    do.call(rbind, lapply(1:nMcmc, function(y) mcmcOut[[y]][[x]]))
  })
  
  return(list(alpha = alphaOut, m = mOut, u = uOut, mcmc = mcmcOut))
}











