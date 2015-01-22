#ifndef __RDFG_h
#define __RDFG_h

#include <Rcpp.h>
#include "rToCpp.h"
#include "PhyDef.h"
#include "DiscreteFactorGraph.h"

using namespace Rcpp;

phy::DFG rToDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors);

class RDFG{
public:
  RDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors);
  double two(){ return 2;}
  int x_;
private:
  phy::DFG dfg;
};

#endif //__RDFG_h