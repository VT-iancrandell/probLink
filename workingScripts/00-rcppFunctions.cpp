#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

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

  // cout << "Source Class: " << sourceClass << endl;
  // cout << "Target Class: "<< targetClass << endl;

  int nLabels = labelsCurrent.size();
  // Here nFeatures doesn't include the indexing columns (the first two)
  int nFeatures =  comparisons.ncol() - 2;


  // Find indices of the records involved in relevant comparisons
  IntegerVector relevantIndices;
  IntegerVector relevantSourceLabels;
  IntegerVector relevantTargetLabels;
  for(int i = 0; i < nLabels; i++){
    //if(i != cppLabelIndex && labelsPropose(i) == targetClass){
    if(labelsPropose(i) == targetClass){
      relevantIndices.push_back(i);
      relevantSourceLabels.push_back(labelsCurrent(i));
      relevantTargetLabels.push_back(labelsPropose(i));
    }
    if(i != cppLabelIndex && labelsCurrent(i) == sourceClass){
    //if(labelsCurrent(i) == sourceClass){
      relevantIndices.push_back(i);
      relevantSourceLabels.push_back(labelsCurrent(i));
      relevantTargetLabels.push_back(labelsPropose(i));
    }
  }
  // cout << "Relevant Indices: "<< relevantIndices << endl;
  // cout << "Proposal labels for SOURCE set: "<< relevantSourceLabels << endl;
  // cout << "Proposal labels for TARGET set: "<< relevantTargetLabels << endl;

  // If there is exactly one relevant index, then the comparison is from a singleton set into another singleton set. This has no effect on the likelihood, so let's just reject this.

  if(relevantIndices.size() == 1){
    //cout << NA_REAL << endl;
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

  //cout << "Row indices from the comparison matrix: "<< relevantComparisonRowIndices << endl;

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
  //cout << "Target posterior:" << targetPosterior << endl;
  //cout << "Source posterior:" << sourcePosterior << endl;
  return targetPosterior - sourcePosterior;
}

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
      // cout << "j =  " << j << endl;
      newLabel = (rand() % n) + 1;
      proposedLabels(j) = newLabel;
      // cout << "Proposed Label =  " << newLabel << endl;
      // cout << "Source labels = " << initialLabels << endl;
      // cout << "Target labels = " << proposedLabels << endl;
      alpha = logProposalRatio(initialLabels, proposedLabels, comparisons, ms, us, priorLinkProb, j + 1);
      // cout << "Alpha = " << alpha << endl;
      if(alpha > log(runif(1)(0))) {
        initialLabels[j] = proposedLabels[j];
        // cout << "Accepted" << endl;
      }else{
        proposedLabels[j] = initialLabels[j];
        // cout << "Rejected" << endl;
      }
    }
    mcmcOut(i,_) = initialLabels;
  }
  return mcmcOut;
}





