#ifndef __Sample_h
#define __Sample_h

#include <iterator>
#include "Definitions.h"
#ifndef BOOST_TEST
#include <Rcpp.h>
#endif

namespace dgRaph {
  template <typename Iterator>
  unsigned discreteSample(Iterator beg, Iterator end){
    #ifndef BOOST_TEST
    // Sample using Runif function
    double u = R::runif(0,1);
    double cumSum = 0;
    double sum = 0;
    Iterator p = beg;

    // Calculate sum
    while(p != end){
      sum += *p++;
    }

    // Find variable
    p = beg;
    while(p != end && cumSum < u*sum){
      cumSum += *p++;
    }

    return std::distance(beg, p) - 1;
    #else
    return 0;
    #endif
  }
}

#endif //__Sample_h
