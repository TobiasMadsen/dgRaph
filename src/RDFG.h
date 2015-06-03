#ifndef __RDFG_h
#define __RDFG_h

#include <Rcpp.h>
#include "rToCpp.h"
#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

phy::DFG rToDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors);

class RDFG{
public:
  // Constructors
  RDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors, IntegerVector potentialMap);

  // Calculate likelihood of of data frame
  double calcLikelihood(IntegerVector observations, LogicalVector observed);
  double calcLogLikelihood(IntegerVector observations, LogicalVector observed);

  // Calculate Expecations
  NumericVector expect(List facScores);

  // Sampling
  IntegerMatrix simulate(int N);

  // Calculate most probable state given partially observed data
  IntegerVector maxProbState(IntegerVector observations, LogicalVector observed);

  // Calculate factor expectation counts
  List facExpCounts(IntegerMatrix observations );

  // Potentials and scores
  List getPotentials();
  void resetPotentials(List facPotentials);
  void resetScores(List facScores);
  
private:
  phy::DFG dfg;

};

#endif //__RDFG_h
