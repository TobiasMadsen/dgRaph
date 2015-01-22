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
  : dfg(rToDFG(varDimensions, facPotentials, facNeighbors)) {}

DataFrame RDFG::makeImportanceSamples(int N, double alpha, List facPotentialsFg){
  //Make conversions
  std::vector<phy::xmatrix_t> facPotFg;
  for(int k = 0; k < facPotentialsFg.size(); ++k){
    facPotFg.push_back( rMatToMat( facPotentialsFg[k] ));
  }

  //Calculate IS distribution
  std::vector<phy::xmatrix_t> facPotIS( dfg.factors.size() );
  for(int f = 0; f < dfg.factors.size(); ++f){
    phy::xmatrix_t const & potNull = dfg.getFactor(f).potential;
    phy::xmatrix_t const & potFg   = facPotFg.at(f);

    phy::xmatrix_t potIS(potNull.size1(), potNull.size2());
    for(int i = 0; i < potIS.size1(); ++i){
      for(int j = 0; j < potIS.size2(); ++j){
	potIS(i,j) = phy::power(potNull(i,j), 1-alpha)*phy::power(potFg(i,j), alpha);
      }
    }
    facPotIS.at(f) = potIS;
  }
  
  phy::DFG dfgFg = dfg;
  phy::DFG dfgIS = dfg;
  dfgFg.resetFactorPotentials( facPotFg );
  dfgIS.resetFactorPotentials( facPotIS );

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
