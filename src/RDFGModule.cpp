#include <Rcpp.h>
#include "RDFG.h"

using namespace Rcpp ;

RCPP_MODULE(phy) {
    class_<RDFG>("RDFG")

      .constructor<IntegerVector, List, List, IntegerVector>()

      .method("simulate", &RDFG::simulate)

      .method("makeImportanceSamples", &RDFG::makeImportanceSamples)

      .method("calculateExpectedScoreIS", &RDFG::calculateExpectedScoreIS)
      
      .method("maxProbState", &RDFG::maxProbState)

      .method("facExpCounts", &RDFG::facExpCounts)

      .method("resetPotentials", &RDFG::resetPotentials)
      
      .method("getPotentials", &RDFG::getPotentials)

      .method("calcLikelihood", &RDFG::calcLikelihood)

      .method("calcLogLikelihood", &RDFG::calcLogLikelihood);

}
