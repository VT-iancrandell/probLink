#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;


// Function to replicate the behavior of x[order(y)]. Also does not modify the arguments in place.

IntegerVector sortBy(IntegerVector x, IntegerVector order){
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){
    return order[i] < order[j];
    });
  return x[idx];
}

//rbindlist version

IntegerMatrix rbindBonds(List pairList){
  
  int nList = pairList.size();
  int nrows = 0;
  for(int i = 0; i < nList; i++){
    nrows += as<IntegerMatrix>(pairList(i)).nrow();
  }
  
  IntegerMatrix out(nrows, 2);
  int listIndex = 0;
  int rowIndex = 0;
  
  for(int i = 0; i < nList; i++){
    IntegerMatrix mat = as<IntegerMatrix>(pairList(listIndex));
    for(int j = 0; j < mat.nrow(); j ++){
      out(rowIndex,_) = mat(j,_);
      rowIndex++;
    }
    listIndex++;
  }
  return out;
}

// Given a set of integers, make all possible pairs of those integers

// [[Rcpp::export]]
NumericMatrix makeCombinations(IntegerVector numbers){
  int n = numbers.size();
  int nc2 = n * (n - 1) / 2;
  NumericMatrix pairs(nc2, 2);
  int row = 0;
  for(int i = 0; i < (n - 1); i++){
    for(int j = (i + 1); j < n; j++){
      pairs(row, 0) = numbers(i);
      pairs(row, 1) = numbers(j);
      row++;
    }
  }
  return pairs;
}

//Given a vector of component labels, make a list of all pairs between rows that have the same label value

// [[Rcpp::export]]
Rcpp::List componentLabelsToBonds(IntegerVector blockingFeature, IntegerVector componentLabels, int nUniqueFeatures){
  int n = componentLabels.size();
  Rcpp::List bondList(nUniqueFeatures);
  
  int startingIndex = 0;
  int endingIndex = 0;
  NumericVector indexSequence;
  int listEntryIndex = 0;
  
  for(int i = 1; i < n; i++){
    if(blockingFeature(i) == blockingFeature(i - 1)){
      endingIndex++;
    }else{
      indexSequence = seq(startingIndex, endingIndex);
      bondList(listEntryIndex) = makeCombinations(componentLabels[indexSequence]);
      startingIndex = endingIndex + 1;
      endingIndex = startingIndex;
      listEntryIndex++;
    }
  }
  indexSequence = seq(startingIndex, endingIndex);
  bondList(listEntryIndex) = makeCombinations(componentLabels[indexSequence]);
  return bondList;
}

// [[Rcpp::export]]
Rcpp::List componentLabelsToBondsCharacter(CharacterVector blockingFeature, IntegerVector componentLabels, int nUniqueFeatures){
  int n = componentLabels.size();
  Rcpp::List bondList(nUniqueFeatures);
  
  int startingIndex = 0;
  int endingIndex = 0;
  NumericVector indexSequence;
  int listEntryIndex = 0;
  
  for(int i = 1; i < n; i++){
    if(blockingFeature(i) == blockingFeature(i - 1)){
      endingIndex++;
    }else{
      indexSequence = seq(startingIndex, endingIndex);
      bondList(listEntryIndex) = makeCombinations(componentLabels[indexSequence]);
      startingIndex = endingIndex + 1;
      endingIndex = startingIndex;
      listEntryIndex++;
    }
  }
  indexSequence = seq(startingIndex, endingIndex);
  bondList(listEntryIndex) = makeCombinations(componentLabels[indexSequence]);
  return bondList;
}

// Finds a bond in a list of bonds (sorted)

// [[Rcpp::export]]
NumericVector findBondInList(int id1, int id2, NumericMatrix bondList){
  int p = bondList.ncol();
  int nrow = bondList.nrow();
  int row = 0;
  
  
  while(id1 != bondList(row, 0) & row < nrow){
    row++;
  }
  while(id2 != bondList(row, 1) & row < nrow){
    row++;
  }
  
  if(row >= bondList.nrow()){
    NumericVector out(p - 2);
    return out;
  }
  
  NumericVector out(p - 2);
  for(int i = 0; i < p - 2; i++){
    out[i] = bondList(row, i+2);
  }
  
  return out;
}

// Functions to draw beta samples for the non-class label parameters

NumericVector sampleBeta(NumericVector a, NumericVector b){
  
  int p = a.size();
  NumericVector out(p);
  
  for(int i = 0; i < p; i++){
    out(i) = rbeta(1, a(i), b(i))(0);
  }
  
  return out;
}

// Takes the valid pairs of rows and makes the 0/1 comparison vectors
// [[Rcpp::export]]
NumericMatrix computeComparisonVectors(IntegerMatrix rowPairs, CharacterMatrix data){
  
  //for m relevant rows, construct the mc2 matrix of comparisons
  int m = rowPairs.nrow();
  int p = data.ncol() - 1;
  NumericMatrix compareLabels(m, p + 2);
  
  //Make vector of comparisons between labels
  for(int i = 0; i < m ; i++){
    
    compareLabels(i, 0) = rowPairs(i, 0);
    compareLabels(i, 1) = rowPairs(i, 1);
    
    for(int k = 0; k < p; k++){
      // Some awkaward indexing here to move around the key column (which comes first). May need to change.
      
      if(CharacterMatrix::is_na(data(rowPairs(i, 0) - 1, k + 1)) || CharacterMatrix::is_na(data(rowPairs(i, 1) - 1, k + 1))){
        compareLabels(i, k + 2) = NA_REAL;
      }else{
        compareLabels(i, k + 2) = data(rowPairs(i, 0) - 1, k + 1) == data(rowPairs(i, 1) - 1, k + 1);
      }
    }
  }
  return compareLabels;
}


// Functions involved in computing the log posterior and the proposal ratio

// [[Rcpp::export]]
float logPosterior(IntegerVector labels, NumericMatrix comparisons, NumericVector ms, NumericVector us, float priorLinkProb){

  int n = labels.size();
  int nc2 = n * (n - 1) / 2;

  //Make vector of comparisons between labels
  int check = 0;
  LogicalVector compareLabels(nc2);
  for(int i = 0; i < n - 1; i++){
    for(int j = i + 1; j < n; j++){
      compareLabels(check) = labels(i) == labels(j);
      check++;
    }
  }
  // Compute rolling sum of comparison by feature terms. If a comparison is missing, disinclude it
  int p = comparisons.ncol();
  float out = 0;
  float sumComponent = 0;
  for(int i = 0; i < nc2; i++){
    NumericVector comparisonRow = comparisons(i,_);
    if(compareLabels(i)){
      for(int j = 0; j < p; j++){
        sumComponent = log(ms(j))*comparisonRow(j) + log(1 - ms(j))*(1 - comparisonRow(j));
        if(!isnan(sumComponent)){
          out += sumComponent;
        }
      }
      out += log(priorLinkProb);
    }else{
      for(int j = 0; j < p; j++){
        sumComponent = log(us(j))*comparisonRow(j) + log(1 - us(j))*(1 - comparisonRow(j));
        if(!isnan(sumComponent)){
          out += sumComponent;
        }
      }
      out += log(1 - priorLinkProb);
    }
  }
  return out;
}

// [[Rcpp::export]]
double logProposalRatio(IntegerVector labelsCurrent, IntegerVector labelsPropose,
                      NumericMatrix comparisons, NumericVector ms, NumericVector us,
                      float priorLinkProb, int changedLabelIndex){

  // Note here that the indexing is being converted from R's 1 indexing to C++'s 0 indexing. May change this later.
  int cppLabelIndex = changedLabelIndex - 1;

  int sourceClass = labelsCurrent(cppLabelIndex);
  int targetClass = labelsPropose(cppLabelIndex);
  int nCompRows = comparisons.nrow();

  int nLabels = labelsCurrent.size();
  
  // Here nFeatures doesn't include the indexing columns (the first two)
  int nFeatures =  comparisons.ncol() - 2;

  // Find indices of the records involved in relevant comparisons
  IntegerVector relevantIndices;
  IntegerVector relevantSourceLabels;
  IntegerVector relevantTargetLabels;
  for(int i = 0; i < nLabels; i++){
    if(labelsPropose(i) == targetClass){
      relevantIndices.push_back(i);
      relevantSourceLabels.push_back(labelsCurrent(i));
      relevantTargetLabels.push_back(labelsPropose(i));
    }
    if(i != cppLabelIndex && labelsCurrent(i) == sourceClass){
      relevantIndices.push_back(i);
      relevantSourceLabels.push_back(labelsCurrent(i));
      relevantTargetLabels.push_back(labelsPropose(i));
    }
  }
  // If there is exactly one relevant index, then the comparison is from a singleton set into another singleton set. This has no effect on the likelihood, so let's just reject this.

  if(relevantIndices.size() == 1){
    //cout << "First na real return" << endl;
    return NA_REAL;
  }

  // Check to see if the needed comparisons exist in the comparison matrix. If they don't, then the blocking assumptions are invalidated and we return NA_REAL

  NumericVector relevantComparisonRowIndices;
  int k = 0;
  bool foundMatch = false;
  int nRelevantIndices = relevantIndices.size();

  for(int i = 0; i < nRelevantIndices - 1; i++){
    for(int j = i + 1; j < nRelevantIndices; j++){
      k = 0;
      while(k < nCompRows){
        if(relevantIndices(i) + 1 == comparisons(k, 0) && relevantIndices(j) + 1 == comparisons(k, 1)){
          relevantComparisonRowIndices.push_back(k);
          foundMatch = true;
          break;
        }
        k++;
      }
      if(!foundMatch){
        return NA_REAL;
      }
      foundMatch = false;
    }
  }


  // Take these row indices and construct the matrix of appropriate comparisons

  NumericMatrix relevantCompMat(relevantComparisonRowIndices.size(), nFeatures);

  for(int i = 0; i < relevantComparisonRowIndices.size(); i++){
    for(int j = 0; j < nFeatures; j++){
      relevantCompMat(i,j) = comparisons(relevantComparisonRowIndices(i), j + 2);
    }
  }

  float targetPosterior;
  float sourcePosterior;
  targetPosterior = logPosterior(relevantTargetLabels, relevantCompMat, ms, us, priorLinkProb);
  sourcePosterior = logPosterior(relevantSourceLabels, relevantCompMat, ms, us, priorLinkProb);
  return targetPosterior - sourcePosterior;
}

// Metropolis algorithm

// [[Rcpp::export]]
NumericMatrix linkageMetropolis(IntegerVector initialLabels, NumericMatrix comparisons, NumericVector ms, NumericVector us, float priorLinkProb, int mcmc, int reportInterval){
  int n = initialLabels.size();
  float alpha;
  int newLabel;

  NumericMatrix mcmcOut(mcmc, n);
  IntegerVector proposedLabels = clone(initialLabels);
  
  for(int i = 0; i < mcmc; i++){
    for(int j = 0; j < n; j++){
      newLabel = (rand() % n) + 1;
      proposedLabels(j) = newLabel;
      alpha = logProposalRatio(initialLabels, proposedLabels, comparisons, ms, us, priorLinkProb, j + 1);
      if(alpha > log(runif(1)(0))) {
        initialLabels[j] = proposedLabels[j];
      }else{
        proposedLabels[j] = initialLabels[j];
      }
    }
    mcmcOut(i,_) = initialLabels;
    if((i + 1) % reportInterval == 0){
      cout << "Beginning iteration " << i << endl;
    }
  }
  return mcmcOut;
}

//
// Functions for extracting the information needed for performing the Gibbs steps
//

// List gibbsStep(IntegerVector labels, NumericMatrix comparisons, NumericVector gammaTotal){
//   
//   int n = labels.size();
//   int nc2 = n * (n - 1) / 2;
//   
//   Rcpp::List bonds = componentLabelsToBonds(sortBy(labels, labels),
//                                             sortBy(seq(1, n), labels),
//                                             unique(labels).size());
//   IntegerMatrix bondMat = rbindBonds(bonds);
//   int nBonds = bondMat.nrow();
//   
//   NumericVector mCounts(comparisons.ncol() - 2);
//   
//   for(int i = 0; i < nBonds; i++){
//     mCounts += findBondInList(bondMat(i, 0), bondMat(i, 1), comparisons);
//   }
//   
//   List out(3);
//   out(0) = rbeta(1, nBonds + 1, nc2 - nBonds + 1)(0);;
//   out(1) = sampleBeta(mCounts + 1, nBonds - mCounts + 1);
//   out(2) = sampleBeta(gammaTotal - mCounts + 1, nc2 - gammaTotal - (nBonds - mCounts) + 1);
//   return out;
//   
// }

// [[Rcpp::export]]  
List gibbsStep(NumericMatrix bondMat, NumericMatrix comparisons, NumericVector gammaTotal){

  int nBonds = bondMat.nrow();

  NumericVector mCounts(comparisons.ncol() - 2);

  for(int i = 0; i < nBonds; i++){
    mCounts += findBondInList(bondMat(i, 0), bondMat(i, 1), comparisons);
  }

  List out(3);
  out(0) = rbeta(1, nBonds + 1, comparisons.nrow() - nBonds + 1)(0);;
  out(1) = sampleBeta(mCounts + 1, nBonds - mCounts + 1);
  out(2) = sampleBeta(gammaTotal - mCounts + 1, comparisons.nrow() - gammaTotal - (nBonds - mCounts) + 1);
  return out;

}







