
# compute the U probabilities from the data

computeU = function(dataCol){
  dataTable = table(dataCol, useNA = 'no')
  N = sum(dataTable)
  sum(dataTable^2 - dataTable) / (N^2 - N)
}

# Use a blocking function to extract the row indices of all pairs of valid comparisons

# extractValidComparisons = function(data, blockingFn){
#   
#   validPairs = list()
#   n = nrow(data)
#   index = 1
#   
#   for(i in 1:(n-1)){
#     for(j in (i + 1):n){
#       if(blockingFn(data[i], data[j])){
#         validPairs[[index]] = c(i, j)
#         index = index + 1
#       }
#     }
#   }
#   
#   return(do.call(rbind, validPairs))
#   
# }

# Use a blocking function to extract the row indices of all pairs of valid comparisons
# If blockingVars == NULL, create the bonds with no blocks

disjunctionBlock = function(blockingVars, data){
  
  if(is.null(blockingVars)) return(makeCombinations(1:nrow(data)))
  
  bonds = vector("list", length(blockingVars))
  
  for(i in 1:length(blockingVars)){
    
    featureIn = data[order(get(blockingVars[i])), list(irow, feature = get(blockingVars[i]))]
    featureInNotSingle = featureIn[,list(irow, .N), by = feature][N > 1]
    bonds[[i]] = do.call(rbind, componentLabelsToBonds(featureInNotSingle$feature, featureInNotSingle$irow, uniqueN(featureInNotSingle$feature)))
    
  }
  
  out = data.table(do.call(rbind, bonds))
  return(as.matrix(unique(out)))
  
}

)