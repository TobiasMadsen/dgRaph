#ifndef __rToCpp_h
#define __rToCpp_h

#include <Rcpp.h>
#include "Definitions.h"
#include "StateMask.h"

using namespace Rcpp;

void dataToStateMasks(IntegerMatrix const & data, unsigned row, dgRaph::stateMaskVec_t & stateMasks);

void dataToStateMasks(IntegerMatrix const & data, List const & dataList, unsigned row, dgRaph::stateMaskVec_t & stateMasks);

void rMatToMat(NumericMatrix const & rmat, dgRaph::matrix_t & mat);

dgRaph::matrix_t rMatToMat(NumericMatrix const & rmat);

std::vector<dgRaph::matrix_t> rFacPotToFacPot(List const & facPot);

List facPotToRFacPot(std::vector<dgRaph::matrix_t> const & facPot);

std::vector<std::vector<unsigned> > rNbsToNbs(List const & rNbs);


#endif  //__rToCpp_h
