# R interface function for the C function computeComparisonVectors

bondsToComparisonVectors = function(bondList, data){
  bonds = as.matrix(bondList)
  data = as.matrix(data)
  out = data.table(computeComparisonVectors(bonds, data))
  colnames(out) = c("id1", "id2", colnames(data)[-1])
  return(out)
}

## Profile the mcmc output

profileMcmc = function(componentMcmc){
  
  mcmcBonds = lapply(1:nrow(componentMcmc), function(x){
    row = sort(componentMcmc[x,])
    out = do.call(rbind, componentLabelsToBondsCharacter(row, as.numeric(names(row)), uniqueN(row)))
    out = out[order(out[,1]),, drop = FALSE]
  })
  
  configurationIndices = numeric(length(mcmcBonds))
  configurationIndices[1] = 1
  
  configurations = vector("list")
  configurationIndex = 1
  configurations[[configurationIndex]] = mcmcBonds[[1]]
  
  for(i in 2:length(mcmcBonds)){
    
    addNewPrototype = TRUE
    
    for(j in 1:configurationIndex){
      if(identical(mcmcBonds[[i]], configurations[[j]])){
        addNewPrototype = FALSE
        break
      }
    }
    if(addNewPrototype){
      configurationIndex = configurationIndex + 1
      configurations[[configurationIndex]] = mcmcBonds[[i]]
      configurationIndices[i] = configurationIndex
    }else{
      configurationIndices[i] = j
    }
  }
  
  configurationFrequencies = data.table(configurationIndex = 1:configurationIndex, frequency = tabulate(configurationIndices)/nrow(componentMcmc))
  return(list(configurations = configurations, whichConfiguration = configurationIndices, configurationFrequencies = configurationFrequencies))
}

# find the most probably configuration from an mcmc profile

findMpc = function(mcmcProfile){
  mpc = do.call(rbind,
                lapply(1:length(mcmcProfile), function(x){
                  mcmcProfile[[x]]$configurations[[which.max(mcmcProfile[[x]]$configurationFrequencies$frequency)]]
                }))
  mpc = data.table(mpc)
  mpcProbs = sapply(1:length(mcmcProfile), function(x) max(mcmcProfile[[x]]$configurationFrequencies$frequency))
  return(list(mpc = mpc, mpcProbs = mpcProbs))
}


