// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include "PhyDef.h"
#include "DiscreteFactorGraph.h"
#include "rToCpp.h"

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

// [[Rcpp::export]]
double mgfDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors){
  std::vector<unsigned> varDim(varDimensions.begin(), varDimensions.end());
  
  std::vector<phy::matrix_t> facPot(facPotentials.size());
  std::vector< std::vector<unsigned> > facNbs(facNeighbors.size()); 
  //Assert same size
  for(int k = 0; k < facPotentials.size(); ++k){
    SEXP lfp = facPotentials[k];
    SEXP lnb = facNeighbors[k];
    
    NumericMatrix mat(lfp);
    IntegerVector nbs(lnb);
    
    phy::matrix_t pot(mat.nrow(), mat.ncol());
    for(int i = 0; i < mat.nrow(); ++i)
      for(int j = 0; j < mat.ncol(); ++j)
        pot(i,j) = mat(i,j);
        
    facPot.at(k) = pot;
    
    std::vector<unsigned> nb;
    for(int i = 0; i < nbs.size(); ++i){
      nb.push_back( nbs(i));
    }
    facNbs.at(k) = nb;
  }
  
  phy::DFG dfg(varDim, facPot, facNbs);
  phy::stateMaskVec_t stateMask( varDim.size() );
  Rcout << dfg.calcNormConst(stateMask) << std::endl;
  return 0;
}

// [[Rcpp::export("PGMExpectCpp")]]
NumericVector PGMExpectCpp(IntegerVector varDimensions, List facPotentials1, List facPotentials2, List facNeighbors){
  //Convert varDimensions to std::vector<unsigned>
  std::vector<unsigned> varDim(varDimensions.begin(), varDimensions.end());
  //Convert both sets of factorPotentials to std::vector<matrix_t>
  std::vector<phy::matrix_t> facPot1;
  std::vector<phy::matrix_t> facPot2;
  for(int k = 0; k < facPotentials1.size(); ++k){
    facPot1.push_back( rMatToMat( facPotentials1[k] ));
    facPot2.push_back( rMatToMat( facPotentials2[k] ));
  }
    
  //Convert factor neighbors to std::vector< std::vector<unsigned> >
  std::vector< std::vector<unsigned> > facNbs = rNbsToNbs( facNeighbors);
  //Calculate fun_a and fun_b
  phy::DFG dfg(varDim, facPot1, facNbs);
  phy::stateMaskVec_t stateMasks(varDimensions.size());
  std::pair<phy::number_t, phy::number_t> res = dfg.calcExpect( facPot1, facPot2, stateMasks);

  NumericVector ret(2);
  ret(0) = res.first;
  ret(1) = res.second;
  return ret;
}
