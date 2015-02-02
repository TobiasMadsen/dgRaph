#include <Rcpp.h>
#include "RDFG.h"

using namespace Rcpp ;

RCPP_MODULE(phy) {
    class_<RDFG>("RDFG")

      .constructor<IntegerVector, List, List>()

      .field("x", &RDFG::x_)

      .method("two", &RDFG::two)

      .method("makeImportanceSamples", &RDFG::makeImportanceSamples)

      .method("calculateExpectedScoreIS", &RDFG::calculateExpectedScoreIS)
      
      .method("maxProbState", &RDFG::maxProbState)

      .method("facExpCounts", &RDFG::facExpCounts)
    ;
}
