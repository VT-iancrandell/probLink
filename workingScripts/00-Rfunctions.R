
# compute the U probabilities from the data

computeU = function(dataCol){
  dataTable = table(dataCol, useNA = 'no')
  N = sum(dataTable)
  sum(dataTable^2 - dataTable) / (N^2 - N)
}
# These functions use precomputed comparison vectors, which require more storage and are potentially faster

logConfigPost = function(labels, comparisons, ms, us, priorLinkProb){
  pairs = combn(labels, 2)
  nc2 = ncol(pairs)
  pairsMatched = apply(pairs, 2, function(x) x[1] == x[2])
  matchProb = numeric(nc2)
  for(j in 1:nc2){
    indicator = unlist(comparisons[j,])
    matchProb[j] = ifelse(pairsMatched[j],
                          sum(log(ms)*indicator + log(1 - ms)*(1 - indicator), na.rm = TRUE) + log(priorLinkProb),
                          sum(log(us)*indicator + log(1 - us)*(1 - indicator), na.rm = TRUE) + log(1 - priorLinkProb))
  }
  sum(matchProb)
}

logProp = function(labelsCurrent, labelsPropose, comparisons, ms, us, priorLinkProb, changedLabelIndex){
  sourceClass = labelsCurrent[changedLabelIndex]
  targetClass = labelsPropose[changedLabelIndex]
  if(targetClass == sourceClass) return(0)

  movedItems = union(sourceClass, targetClass)
  retainedComparisons = which(labelsCurrent %in% c(sourceClass, targetClass))
  relevantComparisons = comparisons[comparisons[,1] %in% retainedComparisons & comparisons[,2] %in% retainedComparisons,-c(1, 2), drop = FALSE]
  if(nrow(relevantComparisons) == 0) return(0)
  relevantNewLabels = labelsPropose[labelsPropose %in% movedItems, drop = FALSE]
  relevantOldLabels = labelsCurrent[labelsCurrent %in% movedItems, drop = FALSE]
  target = logConfigPost(relevantNewLabels, relevantComparisons, ms, us, priorLinkProb)
  source = logConfigPost(relevantOldLabels, relevantComparisons, ms, us, priorLinkProb)
  return(target - source)
}

linkConfigMH = function(labels, comparisons, ms, us, priorLinkProb, mcmc = 1000, reportEvery = 1000){

  n = length(labels)
  mcmcMatrix = matrix(NA, nrow = mcmc, ncol = n)
  labelsPropose = labels

  for(i in 1:mcmc){
    for(j in 1:n){
      labelsPropose[j] = sample(1:n, 1)
      alpha = logProp(labels, labelsPropose, comparisons, ms, us, priorLinkProb, changedLabel = j)
      if(alpha > log(runif(1))) {
        labels[j] = labelsPropose[j]
      }else{
        labelsPropose[j] = labels[j]
      }
    }

    mcmcMatrix[i,] = labels

    if(i %% reportEvery == 0) print(sprintf("Finished iteration %s", i))
  }
  return(mcmcMatrix)
}

linkConfigMHCpp = function(labels, comparisons, ms, us, priorLinkProb, mcmc = 1000, reportEvery = 1000){

  n = length(labels)
  mcmcMatrix = matrix(NA, nrow = mcmc, ncol = n)
  labelsPropose = labels

  for(i in 1:mcmc){
    for(j in 1:n){
      labelsPropose[j] = sample(1:n, 1)
      alpha = logProposalRatio(labels, labelsPropose, comparisons, ms, us, priorLinkProb, changedLabel = j)
      if(is.na(alpha)) alpha = -Inf
      if(alpha > log(runif(1))) {
        labels[j] = labelsPropose[j]
      }else{
        labelsPropose[j] = labels[j]
      }
    }

    mcmcMatrix[i,] = labels

    if(i %% reportEvery == 0) print(sprintf("Finished iteration %s", i))
  }
  return(mcmcMatrix)
}





# These are the same as above but rewritten to compute comparisons on the fly

compareAndValidateBlocking = function(labelsCurrent, labelsPropose, changedLabelIndex, data, checkBlocking = NULL){
  # Check that blocking constraints are satisfied in the target set

  if(!is.null(checkBlocking)){
    targetClass = labelsPropose[changedLabelIndex]
    targetData = data[labelsCurrent == targetClass,]
    check = checkBlocking(data[changedLabelIndex], targetData)
    if(!check){
      return(NULL)
    }
  }

  # Do comparisons given that blocking is satisfied
  targetClass = labelsPropose[changedLabelIndex]
  sourceClass = labelsCurrent[changedLabelIndex]

  relevantData = data[labelsCurrent %in% c(sourceClass, targetClass),]
  nRelevant = nrow(relevantData)
  if(nRelevant == 1) return(NULL)
  ids = t(combn(relevantData$keyby, 2))
  colnames(ids) = c("id1", "id2")
  comparisons = list()
  for(i in 1:nrow(ids)){
    comparisons[[i]] = data[ids[i, 1], -1] == data[ids[i, 2], -1]
  }
  comparisons = do.call(rbind, comparisons)
  return(cbind(ids, comparisons))
}
