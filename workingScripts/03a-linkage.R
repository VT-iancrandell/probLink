source("./src/mcmcLinkReconcile/rcppMCMCDev/02a-recordLinkageMakeBonds.R")


#
# Find non trivial components, create a list with the IDs of the records in each component
#

comps = components(graph_from_edgelist(as.matrix(linkedIds), directed = FALSE))
components = lapply(which(comps$csize > 2), function(x) which(comps$membership == x))
table(sapply(components, length))
length(components)

allComparisonVectors = lapply(components, makeComparisonVectors, demographics, ms = ms, us = us)
priorLinkProb = 1/uniqueN(demographics)


###
### SNOWFALL CODE
###
# Loop over all components and save

nComponents = length(allComparisonVectors)
times = numeric(nComponents)

library(snowfall)
sfInit(parallel=TRUE, cpus=30)
sfExport('allComparisonVectors', 'nComponents', 'times', 'ms', 'us', 'priorLinkProb', 'demographics')
sfSource("./src/mcmcLinkReconcile/00-linkReconcileFunctions.R")

#sfClusterEval(ls())
t0 = Sys.time()

sfSapply(1:300, function(x){
  comparison = allComparisonVectors[[x]]
  nComps = ncol(comparison)
  init = seq_along(comparison$ids)

  #mcmcOut <- linkConfigMH(init, comparison$comparison, ms, us, priorLinkProb, mcmc = 5000, burnin = 500)
  mcmcOut <- linkConfigMH(init, comparison$comparison, ms, us, priorLinkProb, mcmc = 50, burnin = 5)
  colnames(mcmcOut) = as.character(comparison$ids)
  write.csv(mcmcOut, sprintf("./data/dhs_link_char/working/mcmcRecords/wvsMcmcRecords/component%smcmc.csv", x), row.names = FALSE)
  print(sprintf("Ran thing %s.", x))
})

time = Sys.time() - t0



###
### SNOWFALL CODE
###
# Loop over all components and save

nComponents = length(allComparisonVectors)
times = numeric(nComponents)
sfStop
library(snowfall)
sfInit(parallel=TRUE, cpus=30)
sfExport('allComparisonVectors', 'nComponents', 'times', 'ms', 'us', 'priorLinkProb', 'demographics')
sfSource("./src/mcmcLinkReconcile/00-linkReconcileFunctions.R")

#sfClusterEval(ls())
t0 = Sys.time()

sfSapply(1:300, function(x){
  comparison = allComparisonVectors[[x]]
  nComps = ncol(comparison)
  init = seq_along(comparison$ids)

  mcmcOut <- linkConfigMH2(init, comparison$comparison, ms, us, priorLinkProb, mcmc = 5000, burnin = 500)
  colnames(mcmcOut) = as.character(comparison$ids)
  write.csv(mcmcOut, sprintf("./data/dhs_link_char/working/mcmcRecords/wvsMcmcRecords/component%smcmc.csv", x), row.names = FALSE)
  print(sprintf("Ran thing %s.", x))
})

time1 = Sys.time() - t0
time
time1
