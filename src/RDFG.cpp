#include <Rcpp.h>
#include "RDFG.h"
#include "rToCpp.h"
#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

//Helper functions
phy::DFG rToDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors, IntegerVector potentialMap){
    //Make conversions
  std::vector<unsigned> varDim(varDimensions.begin(), varDimensions.end() );
  std::vector<unsigned> potMap(potentialMap.begin(), potentialMap.end() );
  
  std::vector<phy::matrix_t> facPot;
  for(int k = 0; k < facPotentials.size(); ++k){
    facPot.push_back( rMatToMat( facPotentials[k] ));
  }

  std::vector< std::vector<unsigned> > facNbs = rNbsToNbs( facNeighbors);

  return phy::DFG(varDim, facPot, facNbs, potMap);
} 

//Function definition
RDFG::RDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors, IntegerVector potentialMap) 
  : dfg(rToDFG(varDimensions, facPotentials, facNeighbors, potentialMap)) {}

double RDFG::calculateExpectedScoreIS(double alpha, List facPotentialsFg){
  //Make conversions
  std::vector<phy::matrix_t> facPotFg;
  for(int k = 0; k < facPotentialsFg.size(); ++k){
    facPotFg.push_back( rMatToMat( facPotentialsFg[k] ));
  }

  //Calculate IS distribution and score contributions
  std::vector<phy::matrix_t> facPotIS( dfg.factors.size() );
  std::vector<phy::matrix_t> score( dfg.factors.size() );

  for(int f = 0; f < dfg.factors.size(); ++f){
    phy::matrix_t const & potNull = dfg.getFactor(f).getPotential();
    phy::matrix_t const & potFg   = facPotFg.at(f);

    phy::matrix_t potIS(potNull.size1(), potNull.size2());
    phy::matrix_t matScore( potNull.size1(), potNull.size2() );

    for(int i = 0; i < potIS.size1(); ++i){
      for(int j = 0; j < potIS.size2(); ++j){
	potIS(i,j) = phy::power(potNull(i,j), 1-alpha)*phy::power(potFg(i,j), alpha);
	matScore(i,j) = log( potFg(i,j) / potNull(i,j) );
      }
    }
    facPotIS.at(f) = potIS;
    score.at(f) = matScore;
  }

  phy::stateMaskVec_t stateMasks( dfg.variables.size() );
  std::pair<phy::number_t, phy::number_t> res = dfg.calcExpect(facPotIS, score, stateMasks);
  return (res.second / res.first);
}

DataFrame RDFG::makeImportanceSamples(int N, double alpha, List facPotentialsFg){
  //Make conversions
  std::vector<phy::matrix_t> facPotFg;
  for(int k = 0; k < facPotentialsFg.size(); ++k){
    facPotFg.push_back( rMatToMat( facPotentialsFg[k] ));
  }

  //Calculate IS distribution
  std::vector<phy::matrix_t> facPotIS( dfg.factors.size() );
  for(int f = 0; f < dfg.factors.size(); ++f){
    phy::matrix_t const & potNull = dfg.getFactor(f).getPotential();
    phy::matrix_t const & potFg   = facPotFg.at(f);

    phy::matrix_t potIS(potNull.size1(), potNull.size2());
    for(int i = 0; i < potIS.size1(); ++i){
      for(int j = 0; j < potIS.size2(); ++j){
	potIS(i,j) = phy::power(potNull(i,j), 1-alpha)*phy::power(potFg(i,j), alpha);
      }
    }
    facPotIS.at(f) = potIS;
  }
  
  phy::DFG dfgFg = dfg;
  phy::DFG dfgIS = dfg;
  dfgFg.resetPotentials( facPotFg );
  dfgIS.resetPotentials( facPotIS );

  //Calculate normalizing constants for likelihood calculation
  phy::stateMaskVec_t stateMasks( dfg.variables.size() );
  dfg.runSumProduct(stateMasks);
  dfgIS.runSumProduct(stateMasks);
  dfgFg.runSumProduct(stateMasks);
  double normConstNull = dfg.calcNormConst( stateMasks );
  double normConstFg   = dfgFg.calcNormConst( stateMasks );

  dfg.calcFactorMarginals();
  dfg.calcVariableMarginals(stateMasks);
  dfgIS.calcFactorMarginals();
  dfgIS.calcVariableMarginals(stateMasks);

  //Setup
  NumericVector scores(N);
  NumericVector weights(N);
  boost::mt19937 gen(std::time(0));

  const std::vector<phy::vector_t> & varMarNull = dfg.getVariableMarginals();
  const std::vector<phy::matrix_t> & facMarNull = dfg.getFactorMarginals();
  const std::vector<phy::vector_t> & varMarIS   = dfgIS.getVariableMarginals();
  const std::vector<phy::matrix_t> & facMarIS   = dfgIS.getFactorMarginals();

  //Sample
  for(int i = 0; i < N; ++i){
    std::vector<unsigned> sample;
    phy::number_t weight;
    
    dfg.sampleIS(gen,
		 varMarNull,
		 facMarNull,
		 varMarIS,
		 facMarIS,
		 sample,
		 weight);

    phy::number_t scoreNull = dfg.calcFullLikelihood(sample)/normConstNull;
    phy::number_t scoreFg   = dfgFg.calcFullLikelihood(sample)/normConstFg;

    //Output
    scores(i)  = log(scoreFg/scoreNull);
    weights(i) = weight;
  }

  return Rcpp::DataFrame::create(Rcpp::Named("scores")  = scores,
				 Rcpp::Named("weights") = weights );
}

// Preconditions
// observations either empty or length==variables.size()
// observed either empty or length==variables.size()
Rcpp::IntegerVector RDFG::maxProbState(Rcpp::IntegerVector observations, Rcpp::LogicalVector observed){
  //Create empty statemasks
  phy::stateMaskVec_t stateMasks( dfg.variables.size() );
  std::vector<phy::stateMask_t> stateMasksObj;
  stateMasksObj.reserve( dfg.variables.size() );

  //Observed variables
  if(observed.size() == dfg.variables.size()){
    for(int i = 0; i < dfg.variables.size(); ++i){
      if(observed[i]){
	stateMasksObj.push_back( phy::stateMask_t( dfg.getVariable(i).getDimension(), 0 ) );
	// 1-0 R-C++ index conversion
	stateMasksObj.back()( observations[i] - 1 ) = 1;
	stateMasks.at(i) = & stateMasksObj.back();
      }
    }
  }

  //Create return vector
  std::vector<unsigned> maxVarStates( dfg.variables.size() );

  //Calculate most probable state
  dfg.runMaxSum(stateMasks, maxVarStates);
  
  // 0-1 C++-R index conversion
  for(std::vector<unsigned>::iterator it = maxVarStates.begin(); it != maxVarStates.end(); ++it )
    (*it)++;

  return Rcpp::wrap(maxVarStates);
}

Rcpp::List RDFG::facExpCounts(Rcpp::IntegerMatrix observations ){
  if(observations.ncol() != dfg.variables.size() )
    phy::errorAbort("ncol != variables.size");

  dfg.clearCounts();

  //Each row in the matrix is an observation
  for(int i = 0; i < observations.nrow(); ++i){
    //Create statemasks
    phy::stateMaskVec_t stateMasks( dfg.variables.size(), NULL );
    std::vector<phy::stateMask_t> stateMasksObj;
    stateMasksObj.reserve( dfg.variables.size() );

    for(int j = 0; j < observations.ncol(); ++j){
      if( ! IntegerMatrix::is_na(observations(i, j))){
	stateMasksObj.push_back(  phy::stateMask_t( dfg.getVariable(j).getDimension(), 0 ) );
	// 1-0 R-C++ index conversion
	stateMasksObj.back()( observations(i, j) - 1) = 1;
	stateMasks.at(j) = & stateMasksObj.back();
      }
    }

    //calculation
    std::vector<phy::matrix_t> tmpFacMar;
    dfg.initFactorMarginals( tmpFacMar );
    dfg.runSumProduct( stateMasks );
    dfg.calcFactorMarginals( tmpFacMar );
    dfg.submitCounts( tmpFacMar );
  }

  //Convert to list of matrices
  std::vector<phy::matrix_t> potCountsVec;
  dfg.getCounts(potCountsVec);
  Rcpp::List ret(potCountsVec.size());

  for(int p = 0; p < potCountsVec.size(); ++p){
    phy::matrix_t & potCounts = potCountsVec.at(p);
    Rcpp::NumericMatrix rPotCounts(potCounts.size1(), potCounts.size2());
    for(int i = 0; i < potCounts.size1(); ++i)
      for(int j = 0; j < potCounts.size2(); ++j)
	rPotCounts(i,j) = potCounts(i,j);

    ret[p] = rPotCounts;
  }


  return ret;
}

// Accessors
void RDFG::resetPotentials(List facPotentials){
  // Convert to matrix
  std::vector<phy::matrix_t> facPot;
  for(int k = 0; k < facPotentials.size(); ++k){
    facPot.push_back( rMatToMat( facPotentials[k] ));
  }

  // Set factor potentials
  dfg.resetPotentials( facPot);
}

List RDFG::getPotentials(){
  // Loop over factor nodes
  std::vector<phy::matrix_t> ret;
  dfg.getPotentials( ret);
  return facPotToRFacPot( ret);
}


// Preconditions
// observations either empty or length==variables.size()
// observed either empty or length==variables.size()
double RDFG::calcLikelihood(Rcpp::IntegerVector observations, Rcpp::LogicalVector observed){
  // Create empty statemasks
  phy::stateMaskVec_t stateMasks( dfg.variables.size() );
  std::vector<phy::stateMask_t> stateMasksObj;
  stateMasksObj.reserve( dfg.variables.size() );

  // Observed variables
  if(observed.size() == dfg.variables.size()){
    for(int i = 0; i < dfg.variables.size(); ++i){
      if(observed[i]){
	stateMasksObj.push_back( phy::stateMask_t( dfg.getVariable(i).getDimension(), 0 ) );
	// Remember: Take care of 1-0 index-conversion on C++ side
	stateMasksObj.back()( observations[i] - 1) = 1;
	stateMasks.at(i) = & stateMasksObj.back();
      }
    }
  }

  return dfg.calcNormConst(stateMasks);
}
