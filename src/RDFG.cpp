#include <Rcpp.h>
#include "RDFG.h"
#include "rToCpp.h"
#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

//Helper functions
phy::DFG rToDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors){
    //Make conversions
  std::vector<unsigned> varDim(varDimensions.begin(), varDimensions.end() );

  std::vector<phy::xmatrix_t> facPot;
  for(int k = 0; k < facPotentials.size(); ++k){
    facPot.push_back( rMatToMat( facPotentials[k] ));
  }

  std::vector< std::vector<unsigned> > facNbs = rNbsToNbs( facNeighbors);

  return phy::DFG(varDim, facPot, facNbs);
} 

//Function definition
RDFG::RDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors) 
  : dfg(rToDFG(varDimensions, facPotentials, facNeighbors)){}