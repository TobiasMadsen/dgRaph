#include <Rcpp.h>
#include <cmath>
#include <ctime>
#include "RDFG.h"
#include "rToCpp.h"
#include "Definitions.h"
#include "DiscreteFactorGraph.h"
#include "StateMask.h"

using namespace Rcpp;

//Helper functions
dgRaph::DFG rToDFG(IntegerVector const & varDimensions, List const & facPotentials, List const & facNeighbors, IntegerVector const & potentialMap){
    //Make conversions
  std::vector<unsigned> varDim(varDimensions.begin(), varDimensions.end() );
  std::vector<unsigned> potMap(potentialMap.begin(), potentialMap.end() );
  
  std::vector<dgRaph::matrix_t> facPot;
  for(int k = 0; k < facPotentials.size(); ++k){
    facPot.push_back( rMatToMat( facPotentials[k] ));
  }

  std::vector< std::vector<unsigned> > facNbs = rNbsToNbs( facNeighbors);

  return dgRaph::DFG(varDim, facPot, facNbs, potMap);
} 

//Function definitions
RDFG::RDFG(IntegerVector const & varDimensions, List const & facPotentials, List const & facNeighbors, IntegerVector const & potentialMap) 
  : dfg(rToDFG(varDimensions, facPotentials, facNeighbors, potentialMap)) {}

NumericVector RDFG::expect(List const & facScores){
  // Convert to matrix
  std::vector<dgRaph::matrix_t> facScoresVec;
  for(int k = 0; k < facScores.size(); ++k){
    facScoresVec.push_back( rMatToMat( facScores[k] ));
  }

  // Generate statemasks
  dgRaph::stateMaskVec_t stateMasks( dfg.variables.size() );

  // Set factor scores
  dfg.resetScores( facScoresVec);
  std::array<dgRaph::number_t, 2> res = dfg.calcExpect(stateMasks);

  NumericVector ret(2);
  ret(0) = res[0];
  ret(1) = res[1];
  return ret;
}

NumericVector RDFG::gamma(List const & facScores){
  // Convert to matrix
  std::vector<dgRaph::matrix_t> facScoresVec;
  for(int k = 0; k < facScores.size(); ++k){
    facScoresVec.push_back( rMatToMat( facScores[k] ));
  }

  // Generate statemasks
  dgRaph::stateMaskVec_t stateMasks( dfg.variables.size() );

  // Set factor scores
  dfg.resetScores( facScoresVec);
  std::array<dgRaph::number_t, 3> res = dfg.calcGamma(stateMasks);
  
  NumericVector ret(3);
  ret(0) = res[0];
  ret(1) = res[1];
  ret(2) = res[2];
  return ret;
}

IntegerMatrix RDFG::simulate(int N){
  // Calculate marginals
  dgRaph::stateMaskVec_t stateMasks( dfg.variables.size() );
  dfg.runSumProduct(stateMasks);
  dfg.calcVariableMarginals(stateMasks);
  dfg.calcFactorMarginals();
  const std::vector<dgRaph::vector_t> & varMarginals = dfg.getVariableMarginals();
  const std::vector<dgRaph::matrix_t> & facMarginals = dfg.getFactorMarginals();

  // Sample
  IntegerMatrix samples( N, dfg.variables.size() );
  for(int i = 0; i < N; ++i){
    std::vector<unsigned> sample(dfg.variables.size());
    dfg.sample(varMarginals, facMarginals, sample);
    
    // Process
    for(int j = 0; j < sample.size(); ++j)
      samples(i,j) = sample[j] + 1; // C++ to R conversion
  }

  return samples;
}


// Preconditions
// observations either empty or length==variables.size()
// observed either empty or length==variables.size()
IntegerMatrix RDFG::mps(IntegerMatrix const & observations, List const & obsList){
  // Matrix to be returned
  IntegerMatrix ret(observations.nrow(), observations.ncol());

  // Create empty statemasks
  dgRaph::stateMaskVec_t stateMasks( dfg.variables.size() );
  
  for(int i = 0; i < observations.nrow(); ++i){
    // Vector with mps
    std::vector<unsigned> maxVarStates( dfg.variables.size() );

    // Set statemasks
    dataToStateMasks(observations, obsList, i, stateMasks);

    //Calculate most probable state
    dfg.runMaxSum(stateMasks, maxVarStates);

    for(int j = 0; j < observations.ncol(); ++j)
      ret(i, j) = maxVarStates.at(j) + 1;
  }

  return ret;
}

List RDFG::facExpCounts(IntegerMatrix const & observations, List const & obsList ){
  if(observations.ncol() != dfg.variables.size() )
    dgRaph::errorAbort("ncol != variables.size");

  dfg.clearCounts();

  //Create statemasks
  dgRaph::stateMaskVec_t stateMasks( dfg.variables.size());

  //Each row in the matrix is an observation
  for(int i = 0; i < observations.nrow(); ++i){

    dataToStateMasks(observations, obsList, i, stateMasks);

    //calculation
    std::vector<dgRaph::matrix_t> tmpFacMar;
    dfg.initFactorMarginals( tmpFacMar );
    dfg.runSumProduct( stateMasks );
    dfg.calcFactorMarginals( tmpFacMar );
    dfg.submitCounts( tmpFacMar );
  }

  //Convert to list of matrices
  std::vector<dgRaph::matrix_t> potCountsVec;
  dfg.getCounts(potCountsVec);
  List ret(potCountsVec.size());

  for(int p = 0; p < potCountsVec.size(); ++p){
    dgRaph::matrix_t & potCounts = potCountsVec.at(p);
    NumericMatrix rPotCounts(potCounts.size1(), potCounts.size2());
    for(int i = 0; i < potCounts.size1(); ++i)
      for(int j = 0; j < potCounts.size2(); ++j)
	rPotCounts(i,j) = potCounts(i,j);

    ret[p] = rPotCounts;
  }


  return ret;
}

// Accessors
void RDFG::resetPotentials(List const & facPotentials){
  // Convert to matrix
  std::vector<dgRaph::matrix_t> facPot;
  for(int k = 0; k < facPotentials.size(); ++k){
    facPot.push_back( rMatToMat( facPotentials[k] ));
  }

  // Set factor potentials
  dfg.resetPotentials( facPot);
}

void RDFG::resetScores(List const & facScores){
  // Convert to matrix
  std::vector<dgRaph::matrix_t> facScoresVec;
  for(int k = 0; k < facScores.size(); ++k){
    facScoresVec.push_back( rMatToMat( facScores[k] ));
  }

  // Set factor potentials
  dfg.resetPotentials( facScoresVec);
}

List RDFG::getPotentials(){
  // Loop over factor nodes
  std::vector<dgRaph::matrix_t> ret;
  dfg.getPotentials( ret);
  return facPotToRFacPot( ret);
}

NumericVector RDFG::calcLikelihood(IntegerMatrix const & observations, List const & obsList){
  return exp( calcLogLikelihood(observations, obsList));
}

NumericVector RDFG::calcLogLikelihood(IntegerMatrix const & observations, List const & obsList){
  // Vector to be returned
  NumericVector ret( observations.nrow() );
  
  // Create empty statemasks
  dgRaph::stateMaskVec_t stateMasks;

  // Observed variables
  for(int i = 0; i < observations.nrow(); ++i){
    dataToStateMasks(observations, obsList, i, stateMasks);
    ret(i) = dfg.calcLogNormConst(stateMasks);
  }

  return ret;
}
