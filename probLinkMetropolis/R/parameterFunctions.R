#
# This file contains fuctions used to compute and set initial parameter configurations
#

# compute the U probabilities from the data

computeU = function(dataCol){
  dataTable = table(dataCol, useNA = 'no')
  N = sum(dataTable)
  sum(dataTable^2 - dataTable) / (N^2 - N)
}

# An easy function for default parameters. M is a series of .75, U is computed assuming no links, the prior is 1/nrow(data), and runs are set to 1000

defaultParameters = function(data){
  out = list(ms = rep(.75, ncol(data) - 1),
       us = sapply(data[,-1], computeU),
       priorLinkProb = 1/nrow(data),
       nMcmc = 1000)
  return(out)
}