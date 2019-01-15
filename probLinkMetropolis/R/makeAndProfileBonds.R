
# Bond function based on adist (levenstein) for fuzzy matching

makeFuzzyLink = function(fuzzyCols, irow, dmax = 1){
  
  # This if statement links everything if no fuzzy column is provided
  
  if(is.null(fuzzyCols)) {
    fuzzyAdj = matrix(TRUE, length(irow), length(irow))
  }else{
    fuzzyAdj = adist(as.character(fuzzyCols)) <= dmax
  }
  
  g <- graph_from_adjacency_matrix(fuzzyAdj, "upper", diag = F)
  out = as.data.table(matrix(irow[get.edgelist(g)], ncol = 2))
  
  return(out)
}

# makeBonds takes a vector of column names for exact matches and a single character for fuzzy matches. The bonds are tracked by an indicator column, internally labelled irow. 

makeBonds = function(dt, exactCols = NULL, fuzzyCols = NULL, dmax = 1){
  
  # Add indicator for row if it doesn't exist.
  if(!("irow" %in% colnames(dt))) dt$irow = (1:nrow(dt))
  
  # Giving NULL for both the exact columns and fuzzycols links everything
  
  if(is.null(exactCols) & is.null(fuzzyCols)){
    bondsDT = na.omit(dt[, 
                         makeFuzzyLink(NULL, irow)])
  }else if(!is.null(exactCols) & is.null(fuzzyCols)){
    bondsDT = na.omit(dt[, 
                         makeFuzzyLink(NULL, irow), 
                         keyby = c(exactCols)],
                      cols = exactCols)
  }else{
    bondsDT = na.omit(dt[, 
                         makeFuzzyLink(get(fuzzyCols), irow), 
                         keyby = c(exactCols)],
                      cols = exactCols)
  }
  return(bondsDT[,.(V1, V2)])
  
}

# Listified version of the above

makeBondsFromList = function(records, exactRules = NULL, fuzzyRules = NULL, dmax = 1){
  
  if(is.null(exactRules) & is.null(fuzzyRules)){
    out = makeBonds(records)
  }else{
    allBonds = mapply(makeBonds, exactCols = exactRules, fuzzyCols = fuzzyRules, MoreArgs = list(dt = records, dmax = dmax), SIMPLIFY = FALSE)
    out = unique(do.call('rbind', allBonds))
  }
  return(out)
}

# disjunctionBlock = function(blockingVars, data){
# 
#   if(is.null(blockingVars)) return(makeCombinations(1:nrow(data)))
# 
#   bonds = vector("list", length(blockingVars))
# 
#   for(i in 1:length(blockingVars)){
# 
#     featureIn = data[order(get(blockingVars[i])), list(irow, feature = get(blockingVars[i]))]
#     featureInNotSingle = featureIn[,list(irow, .N), by = feature][N > 1]
#     if(nrow(featureInNotSingle) == 0){
#       bonds[[i]] = NULL
#     }else{
#       bonds[[i]] = do.call(rbind, componentLabelsToBonds(featureInNotSingle$feature, featureInNotSingle$irow, uniqueN(featureInNotSingle$feature)))
#     }
#   }
#   out = data.table(do.call(rbind, bonds))
#   return(as.matrix(unique(out)))
# }



#
# Profiling code
#

# This function treats a bond list as an edge list, makes a graph, and gives its component statistics

componentSizesFromBonds = function(bondList){
  #if(nrow(bondList == 0)) return(NA)
  
  g = graph_from_edgelist(as.matrix(bondList), directed = FALSE)
  return(components(g))
}

# This computes the precision and recall for a blocking scheme given a true labelling.

profileBonds = function(trueLabels, blockingBonds, data){
  
  trueBonds = makeBonds(data.table(trueLabels), "trueLabels")
  setkey(trueBonds, V1, V2)
  setkey(blockingBonds, V1, V2)
  
  # trueBonds$foundByMcmc = FALSE
  # blockingBonds$isTrueBond = FALSE
  
  truePositives = copy(trueBonds)[blockingBonds, nomatch = 0]
  falseNegatives = copy(trueBonds)[!blockingBonds]
  falsePositives = blockingBonds[!trueBonds]
  
  tpCount = nrow(truePositives)
  fpCount = nrow(falsePositives)
  
  precision = tpCount / (tpCount + fpCount)
  recall = tpCount / nrow(trueBonds)
  
  colnames(truePositives) = c("id1", "id2")
  colnames(falsePositives) = c("id1", "id2")
  colnames(falseNegatives) = c("id1", "id2")
  
  # Get records that were missed
  
  missed = data[irow %in% unlist(falseNegatives)]
  
  missed = data.table(missed, trueLabel = trueLabels[missed$irow], key = 'trueLabel')
  
  out = list(truePositives = truePositives, falsePositives = falsePositives, falseNegatives = falseNegatives, precision = precision, recall = recall, missedLinks = missed)
  return(out)
}











