#ifndef __rToCpp_h
#define __rToCpp_h

#include <Rcpp.h>
#include "PhyDef.h"

using namespace Rcpp;

void rMatToMat(NumericMatrix const & rmat, phy::xmatrix_t & mat);

phy::xmatrix_t rMatToMat(NumericMatrix const & rmat);

std::vector<std::vector<unsigned> > rNbsToNbs(List const & rNbs);

#endif  //__rToCpp_h
