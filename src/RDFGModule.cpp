#include <Rcpp.h>
#include "RDFG.h"

using namespace Rcpp ;

RCPP_MODULE(phy) {
    class_<RDFG>("RDFG")

      .constructor<IntegerVector, List, List, IntegerVector>()

      .method("makeImportanceSamples", &RDFG::makeImportanceSamples)

      .method("calculateExpectedScoreIS", &RDFG::calculateExpectedScoreIS)
      
      .method("maxProbState", &RDFG::maxProbState)

      .method("facExpCounts", &RDFG::facExpCounts)

      .method("resetFactorPotentials", &RDFG::resetFactorPotentials)
      
      .method("getFactorPotentials", &RDFG::getFactorPotentials)

      .method("calcLikelihood", &RDFG::calcLikelihood);
    ;
}
