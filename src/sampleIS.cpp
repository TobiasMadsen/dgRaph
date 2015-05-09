#include <Rcpp.h>
#include "rToCpp.h"
#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sampleISCpp(int N, double alpha, IntegerVector varDimensions, List facPotentialsNull, List facNeighbors, List facPotentialsFg) {
  //Make conversions
  std::vector<unsigned> varDim(varDimensions.begin(), varDimensions.end() );

  std::vector<phy::matrix_t> facPotNull;
  std::vector<phy::matrix_t> facPotFg;
  for(int k = 0; k < facPotentialsNull.size(); ++k){
    facPotNull.push_back( rMatToMat( facPotentialsNull[k] ));
    facPotFg.push_back( rMatToMat( facPotentialsFg[k] ));
  }

  std::vector< std::vector<unsigned> > facNbs = rNbsToNbs( facNeighbors);

  Rcout << "Debug1" << std::endl;
  phy::DFG dfgNull(varDim, facPotNull, facNbs);
  Rcout << "Debug2" << std::endl;
  phy::DFG dfgFg(varDim, facPotFg, facNbs);

  //Calculate IS distribution
  std::vector<phy::matrix_t> facPotIS( facPotNull.size() );
  for(int f = 0; f < facPotNull.size(); ++f){
    phy::matrix_t const & potNull = facPotNull.at(f);
    phy::matrix_t const & potFg   = facPotFg.at(f);

    phy::matrix_t potIS(potNull.size1(), potNull.size2());
    for(int i = 0; i < potIS.size1(); ++i){
      for(int j = 0; j < potIS.size2(); ++j){
	potIS(i,j) = phy::power(potNull(i,j), 1-alpha)*phy::power(potFg(i,j), alpha);
      }
    }
    facPotIS.at(f) = potIS;
  }
  phy::DFG dfgIS(varDim, facPotIS, facNbs);

  
  //Calculate normalizing constants for likelihood calculation
  phy::stateMaskVec_t stateMasks( varDim.size() );
  dfgNull.runSumProduct(stateMasks);
  dfgIS.runSumProduct(stateMasks);
  dfgFg.runSumProduct(stateMasks);
  double normConstNull = dfgNull.calcNormConst( stateMasks );
  double normConstFg   = dfgFg.calcNormConst( stateMasks );

  dfgNull.calcFactorMarginals();
  dfgNull.calcVariableMarginals(stateMasks);
  dfgIS.calcFactorMarginals();
  dfgIS.calcVariableMarginals(stateMasks);

  //Setup
  NumericVector scores(N);
  NumericVector weights(N);
  boost::mt19937 gen(std::time(0));

  const std::vector<phy::vector_t> & varMarNull = dfgNull.getVariableMarginals();
  const std::vector<phy::matrix_t> & facMarNull = dfgNull.getFactorMarginals();
  const std::vector<phy::vector_t> & varMarIS   = dfgIS.getVariableMarginals();
  const std::vector<phy::matrix_t> & facMarIS   = dfgIS.getFactorMarginals();

  //Sample
  for(int i = 0; i < N; ++i){
    std::vector<unsigned> sample;
    phy::number_t weight;
    
    dfgNull.sampleIS(gen,
		     varMarNull,
		     facMarNull,
		     varMarIS,
		     facMarIS,
		     sample,
		     weight);

    phy::number_t scoreNull = dfgNull.calcFullLikelihood(sample)/normConstNull;
    phy::number_t scoreFg   = dfgFg.calcFullLikelihood(sample)/normConstFg;

    //Output
    scores(i)  = log(scoreFg/scoreNull);
    weights(i) = weight;
  }

  return Rcpp::DataFrame::create(Rcpp::Named("scores")  = scores,
				 Rcpp::Named("weights") = weights );
}
