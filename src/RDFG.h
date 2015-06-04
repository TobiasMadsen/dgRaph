#ifndef __RDFG_h
#define __RDFG_h

#include <Rcpp.h>
#include "rToCpp.h"
#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

phy::DFG rToDFG(IntegerVector const & varDimensions, List const & facPotentials, List const & facNeighbors);

class RDFG{
public:
  // Constructors
  RDFG(IntegerVector const & varDimensions, List const & facPotentials, List const & facNeighbors, IntegerVector const & potentialMap);

  // Calculate likelihood of of data frame
  NumericVector calcLikelihood(IntegerMatrix const & observations);
  NumericVector calcLogLikelihood(IntegerMatrix const & observations);

  // Calculate Expecations
  NumericVector expect(List const & facScores);

  // Sampling
  IntegerMatrix simulate(int N);

  // Calculate most probable state given partially observed data
  IntegerMatrix mps(IntegerMatrix const & observations);

  // Calculate factor expectation counts
  List facExpCounts(IntegerMatrix const & observations );

  // Potentials and scores
  List getPotentials();
  void resetPotentials(List const & facPotentials);
  void resetScores(List const & facScores);
  
private:
  phy::DFG dfg;

};

#endif //__RDFG_h
