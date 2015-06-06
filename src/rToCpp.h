#ifndef __rToCpp_h
#define __rToCpp_h

#include <Rcpp.h>
#include "PhyDef.h"
#include "StateMask.h"

using namespace Rcpp;

void dataToStateMasks(IntegerMatrix const & data, unsigned row, phy::stateMaskVec_t & stateMasks);

void dataToStateMasks(IntegerMatrix const & data, List const & dataList, unsigned row, phy::stateMaskVec_t & stateMasks);

void rMatToMat(NumericMatrix const & rmat, phy::matrix_t & mat);

phy::matrix_t rMatToMat(NumericMatrix const & rmat);

std::vector<phy::matrix_t> rFacPotToFacPot(List const & facPot);

List facPotToRFacPot(std::vector<phy::matrix_t> const & facPot);

std::vector<std::vector<unsigned> > rNbsToNbs(List const & rNbs);


#endif  //__rToCpp_h
