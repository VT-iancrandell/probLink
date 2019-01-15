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
  
  partitionedComparisonVectors = lapply(1:nComps, function(x) {
    
    initialLabels = 1:length(noSingletons[[x]])
    activeComp = comparisonVectors[id1 %in% noSingletons[[x]]]
    relabeling = cbind(match(activeComp$id1, noSingletons[[x]]), match(activeComp$id2, noSingletons[[x]]))
    out = linkageMetropolis(initialLabels, as.matrix(cbind(relabeling, activeComp[,-c(1, 2)])), ms, us, priorLinkProb, mcmc = nMcmc, reportInterval = reportInterval)
    colnames(out) = noSingletons[[x]]
    apply(out, 2, function(y) paste0(y, ".", x))
  })
  
}







# Look at probable configurations and graps

pairwiseEquality = function(labels){
  comparisons = t(combn(labels, 2))
  pairsMatched = apply(comparisons, 1, function(x) x[1] == x[2])
  return(pairsMatched)
}

# Makes an equivalence class matrix. An equivalence class is defined by its pairwise equalities, and so is invariant to label permutations.
# This needs to be redone to ignore non links, so rather than a class being a logical(choose(n, 2)) vector, it's just an 'edgelist'.

makeEcm = function(mcmcOut){
  pairs = t(apply(mcmcOut, 1, pairwiseEquality))
  if(nrow(pairs) == 1){
    equivClassMatrix = data.table(t(pairs))
  }else{
    equivClassMatrix = data.table(pairs)
  }
  colnames(equivClassMatrix) = apply(t(combn(colnames(mcmcOut), 2)), 1, paste0, collapse = 'to')
  return(equivClassMatrix[,.N, by = c(colnames(equivClassMatrix))][order(-N)])
}

ecmTabToGraph = function(ecmTab){
  ncols = ncol(ecmTab)
  allPairs = do.call(rbind, strsplit(colnames(ecmTab)[-ncols], "to"))
  allVert = unique(c(allPairs))
  outGraphs = list()
  for(i in 1:nrow(ecmTab)){
    outGraphs[[i]] = graph_from_edgelist(allPairs[unlist(ecmTab[i,-ncols, with = FALSE]), ,drop = FALSE], directed = FALSE)
    singletons = setdiff(allVert, as.character(names(V(outGraphs[[i]]))))
    outGraphs[[i]] = add_vertices(outGraphs[[i]], length(singletons), name = singletons)
  }
  return(outGraphs)
}

# R interface function for the C function computeComparisonVectors

bondsToComparisonVectors = function(bondList, data){
  bonds = as.matrix(bondList)
  data = as.matrix(data)
  out = data.table(computeComparisonVectors(bonds, data))
  colnames(out) = c("id1", "id2", colnames(data)[-1])
  return(out)
}

## Profile bonds code

profileMcmcBonds = function(trueLabels, mcmcLabels){
  
  nRuns = nrow(mcmcLabels)
  precision = numeric(nRuns)
  recall = numeric(nRuns)
  foundBondsList = vector("list", nRuns)
  falsePositiveList = vector("list", nRuns)
  
  trueBonds = data.table(t(apply(disjunctionBlock('trueLabels', data.table(irow = 1:length(trueLabels), trueLabels = trueLabels)), 1, sort)))
  setkey(trueBonds, V1, V2)
  
  for(i in 1:nRuns){
    mcmcBonds = data.table(t(apply(disjunctionBlock('mcmcLabels', data.table(irow = 1:length(trueLabels), mcmcLabels = mcmcLabels[i,])), 1, sort)))
    if(nrow(mcmcBonds) == 0){
      next
    }
    setkey(mcmcBonds, V1, V2)
    
    # trueBonds[, foundByMcmc := FALSE][mcmcBonds, foundByMcmc := TRUE]
    # mcmcBonds[, isTrueBond := FALSE][trueBonds, isTrueBond := TRUE]
    
    trueBonds$foundByMcmc = FALSE
    mcmcBonds$isTrueBond = FALSE
    trueBonds[mcmcBonds, foundByMcmc := TRUE]
    mcmcBonds[trueBonds, isTrueBond := TRUE]
    
    precision[i] = sum(mcmcBonds$isTrueBond)/nrow(mcmcBonds)
    recall[i] = sum(trueBonds$foundByMcmc)/nrow(trueBonds)
    foundBondsList[[i]] = trueBonds$foundByMcmc
    if(mean(mcmcBonds$isTrueBond) < 1) falsePositiveList[[i]] = data.table(mcmcBonds[isTrueBond == FALSE, list(V1, V2)], runIndex = i)
  }
  
  foundBonds = data.table(trueBonds[,list(V1, V2)], do.call('cbind', foundBondsList))
  colnames(foundBonds) = c("id1", "id2", paste0("Run", 1:nRuns))
  falsePositives = do.call('rbind', falsePositiveList)
  colnames(falsePositives) = c("id1", "id2", 'runIndex')
  
  out = list(recallMatrix = foundBonds, falsePositives = falsePositives, precision = precision, recall = recall)
  return(out)
}


# find most probable configuration and its probability

findMPC = function(mcmcOut, includePairwiseProbabilities = FALSE, report = TRUE){
  
  nComps = length(mcmcOut)
  probsOfMPC = numeric(nComps)
  mpcs = vector('list', nComps)
  pairwiseProbs = vector('list', nComps)
  
  for(i in 1:nComps){
    ecm = makeEcm(mcmcOut[[i]])
    ncols = ncol(ecm)
    mpcs[[i]] = do.call('rbind',strsplit(colnames(ecm[1])[-ncols], "to"))[unlist(ecm[1,-ncols, with = FALSE]),, drop = FALSE]
    probsOfMPC[i] = as.numeric(ecm[1, ..ncols]/sum(ecm[,..ncols]))
    if(includePairwiseProbabilities){
      counts = apply(as.matrix(ecm)[,-ncols, drop = FALSE], 2, function(x) x * as.matrix(ecm)[,ncols, drop = FALSE])
      
      if(is.null(dim(counts))){
        pairwiseProbs[[i]] = data.table(do.call('rbind',strsplit(colnames(ecm)[-ncols], "to")), pLink = counts/sum(ecm[,..ncols]))
      }else{
        counts = colSums(apply(as.matrix(ecm)[,-ncols, drop = FALSE], 2, function(x) x * as.matrix(ecm)[,ncols, drop = FALSE]))
        pairwiseProbs[[i]] = data.table(do.call('rbind',strsplit(colnames(ecm)[-ncols], "to")), pLink = counts/sum(ecm[,..ncols]))
      }
    }
    if(report) print(paste("Finished component", i))
  }
  
  mpcSet = data.table(do.call(rbind, mpcs))
  mpcSet[,V1 := as.numeric(V1)]
  mpcSet[,V2 := as.numeric(V2)]
  mpcProbs = data.table(componentId = 1:nComps, probOfMPC = probsOfMPC)
  
  if(includePairwiseProbabilities){
    pairwiseProbs = do.call('rbind', pairwiseProbs)
    colnames(pairwiseProbs) = c("id1", "id2", "pLink")
    pairwiseProbs = pairwiseProbs[pLink > 0]
    return(list(componentID = mpcSet, mpcProbs = mpcProbs, pairwiseProbs = pairwiseProbs))
  }else{
    return(list(componentID = mpcSet, mpcProbs = mpcProbs))
  }
}




