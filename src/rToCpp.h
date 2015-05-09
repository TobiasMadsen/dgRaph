#ifndef __rToCpp_h
#define __rToCpp_h

#include <Rcpp.h>
#include "PhyDef.h"

using namespace Rcpp;

void rMatToMat(NumericMatrix const & rmat, phy::matrix_t & mat);

phy::matrix_t rMatToMat(NumericMatrix const & rmat);

std::vector<phy::matrix_t> rFacPotToFacPot(List const & facPot);

List facPotToRFacPot(std::vector<phy::matrix_t> const & facPot);

std::vector<std::vector<unsigned> > rNbsToNbs(List const & rNbs);


#endif  //__rToCpp_h
