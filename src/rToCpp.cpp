#include <Rcpp.h>
#include "rToCpp.h"

using namespace Rcpp;

void dataToStateMasks(IntegerMatrix const & data, unsigned row, phy::stateMaskVec_t & stateMasks){
  stateMasks.resize( data.ncol() );
  for(int col = 0; col < data.ncol(); ++col){
    stateMasks.at(col).reset();
    if( ! IntegerMatrix::is_na( data(row, col)) ){
      stateMasks.at(col) = phy::stateMaskPtr_t( new phy::StateMaskObserved( data(row, col) - 1));
    }
  }
}

void dataToStateMasks(IntegerMatrix const & data, List const & dataList, unsigned row, phy::stateMaskVec_t & stateMasks){
  stateMasks.resize( data.ncol() );
  for(int col = 0; col < data.ncol(); ++col){
    stateMasks.at(col).reset();

    if( col < dataList.size() && ! Rf_isNull( dataList[col]) ){
      // Use dataList if possible
      NumericVector x = ((List) dataList[col])[row];
      phy::vector_t posterior(x.size());
      for(int i = 0; i < x.size(); ++i)
	posterior(i) = x(i);
      stateMasks.at(col) = phy::stateMaskPtr_t( new phy::StateMaskPosterior( posterior) );
      continue;
    }

    if( ! IntegerMatrix::is_na( data(row, col)) ){
      // Use observed variable
      stateMasks.at(col) = phy::stateMaskPtr_t( new phy::StateMaskObserved( data(row, col) - 1));
    }

  }
}

void rMatToMat(NumericMatrix const & rmat, phy::matrix_t & mat){
  mat.resize( rmat.nrow(), rmat.ncol());
  for(int i = 0; i < rmat.nrow(); ++i)
    for(int j = 0; j < rmat.ncol(); ++j)
      mat(i,j) = rmat(i,j);
}

phy::matrix_t rMatToMat(NumericMatrix const & rmat){
  phy::matrix_t ret;
  rMatToMat(rmat, ret);
  return ret;
}

std::vector<phy::matrix_t> rFacPotToFacPot(List const & facPot){
  std::vector<phy::matrix_t> ret;
  for(int k = 0; k < facPot.size(); ++k){
    ret.push_back( rMatToMat( facPot[k] ));
  }

  return ret;
}

List facPotToRFacPot(std::vector<phy::matrix_t> const & facPot){
  List ret(facPot.size());
  int fcount = 0;
  for(std::vector<phy::matrix_t>::const_iterator it = facPot.begin(); it != facPot.end(); ++it){
    NumericMatrix m( it->size1(), it->size2() );
    for(int i = 0; i < it->size1(); ++i){
      for(int j = 0; j < it->size2(); ++j){
	m(i,j) = (*it)(i,j);
      }
    }
    ret[fcount] = m;
    fcount++;
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
