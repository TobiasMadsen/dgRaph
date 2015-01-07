// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

// [[Rcpp::export]]
double evalDFG(NumericVector x) {
   std::vector<unsigned> varDim;
   varDim.push_back(x.size());
   
   std::vector<phy::matrix_t> facPot(1);
   facPot.at(0) = phy::matrix_t(1,x.size());
   for(int i = 0; i < x.size(); ++i)
      facPot.at(0)(0,i) = x(i);
      
   std::vector< std::vector<unsigned> > facNbs(1);
   facNbs.at(0).push_back(0);
   
   phy::DFG dfg(varDim, facPot, facNbs);
   
   phy::stateMaskVec_t stateMasks(1);
   return dfg.calcNormConst(stateMasks);
}

