#include <Rcpp.h>
#include "RDFG.h"

using namespace Rcpp ;

RCPP_MODULE(phy) {
    class_<RDFG>("RDFG")

    .constructor<IntegerVector, List, List>()

    .field("x", &RDFG::x_)

    .method("two", &RDFG::two)
    ;
}
