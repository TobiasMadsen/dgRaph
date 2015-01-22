#include <Rcpp.h>
#include "rToCpp.h"

using namespace Rcpp;

void rMatToMat(NumericMatrix const & rmat, phy::xmatrix_t & mat){
  mat.resize( rmat.nrow(), rmat.ncol());
  for(int i = 0; i < rmat.nrow(); ++i)
    for(int j = 0; j < rmat.ncol(); ++j)
      mat(i,j) = rmat(i,j);
}

phy::xmatrix_t rMatToMat(NumericMatrix const & rmat){
  phy::xmatrix_t ret;
  rMatToMat(rmat, ret);
  return ret;
}

std::vector<phy::xmatrix_t> rFacPotToFacPot(List const & facPot){
  std::vector<phy::xmatrix_t> ret;
  for(int k = 0; k < facPot.size(); ++k){
    ret.push_back( rMatToMat( facPot[k] ));
  }

  return ret;
}

std::vector<std::vector<unsigned> > rNbsToNbs(List const & rNbs){
  std::vector< std::vector<unsigned> > ret(rNbs.size());

  for(int i = 0; i < rNbs.size(); ++i){

    SEXP sexpNbs = rNbs[i];
    IntegerVector IVnbs(sexpNbs);
    std::vector<unsigned> nbs;
    for(int j = 0; j < IVnbs.size(); ++j)
      nbs.push_back( IVnbs(j) - 1 );
    ret.at(i) = nbs;
  }
  return ret;
}
