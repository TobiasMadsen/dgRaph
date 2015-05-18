#define BOOST_TEST_DYN_LINK

#include <cmath>

// tested code ( the defines allow access to private data
// consider if test should only include interface not implementation
#define private public
#define protected public
#include "../../src/DiscreteFactorGraph.h"
#undef private
#undef protected

// bost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using boost::unit_test::test_suite;
using namespace phy;

// Constants
#define EPS 1e-10

/**************************************************
 Helper Functions
 **************************************************/

DFG mkSimpleFactorGraph()
{
  // strategy: define simple factor graph and check properties
  // factor graph is a phylo tree with one root and two leaves (and one initial distribution factor).
  vector<unsigned> varDimensions;
  vector<matrix_t> facPotentials;
  vector<vector<unsigned> > varNeighbors(3);
  vector<vector<unsigned> > facNeighbors(3);
  
  // three variables
  unsigned dim = 4;
  varDimensions.push_back(dim);
  varDimensions.push_back(dim);
  varDimensions.push_back(dim);
  
  // 1D-potential
  matrix_t pot1D(1, dim);
  for (unsigned i = 0; i < dim; i++)
    pot1D(0, i) = 0.25;
  
  // 2D-potential
  matrix_t pot2D(dim, dim);
  for(unsigned i = 0; i < dim; ++i){
    for(unsigned j = 0; j < dim; ++j){
      pot2D(i,j) = 0.1;
      if( i == j)
	pot2D(i,j) = 0.7;
    }
  }

  // three factors
  facPotentials.push_back(pot2D);  // pot2D = P(x,y) = exp(Qt)
  facPotentials.push_back(pot2D);
  facPotentials.push_back(pot1D);    // pot1D = equiFreq = [0.25, 0.25, 0.25, 0.25]
  
  // setting up the graph (tree) structure
  // this is done purely in terms of factor links

  // fac 0
  facNeighbors[0].push_back(0);  // neighbor order defines variable order for factor potentials!
  facNeighbors[0].push_back(1);

  // fac 1
  facNeighbors[1].push_back(0);  // neighbor order defines variable order for factor potentials!
  facNeighbors[1].push_back(2);

  // fac 2
  facNeighbors[2].push_back(0);
  
  return DFG(varDimensions, facPotentials, facNeighbors);
}

/**************************************************
 Tests
 **************************************************/

BOOST_AUTO_TEST_CASE(Test_1) 
{
  BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(DFGNode_1) 
{
  // test variable properties
  unsigned dim = 4;
  DFGNode varNode(dim);
  
  BOOST_CHECK( varNode.getDimension() == dim );
  BOOST_CHECK( varNode.isFactor() == false );
  BOOST_CHECK( varNode.getPotential().size1() == 0);
  BOOST_CHECK( varNode.getPotential().size2() == 0);
}


BOOST_AUTO_TEST_CASE(DFGNode_2) 
{
  // test 1D-factor properties
  unsigned size = 4;
  matrix_t potential(1, size);

  Potential pot(potential);
  DFGNode facNode(&pot);
  
  BOOST_CHECK( facNode.getDimension() == 1 );
  BOOST_CHECK( facNode.isFactor() == true );
  BOOST_CHECK( facNode.getPotential().size1() == 1);
  BOOST_CHECK( facNode.getPotential().size2() == size);
}


BOOST_AUTO_TEST_CASE(DFGNode_3) 
{
  // test 2D-factor properties
  unsigned size = 4;
  matrix_t potential(size, size,0);

  Potential pot(potential);
  DFGNode facNode(&pot);
  
  BOOST_CHECK( facNode.getDimension() == 2 );
  BOOST_CHECK( facNode.isFactor() == true );
  BOOST_CHECK( facNode.getPotential().size1() == size);
  BOOST_CHECK( facNode.getPotential().size2() == size);
}

BOOST_AUTO_TEST_CASE(DFG_1) 
{
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)

  // check that links between nodes are defined both ways. 
  for (unsigned i = 0; i < fg.neighbors.size(); i++)
    for (unsigned j = 0; j < fg.neighbors[i].size(); j++) {
      unsigned other = fg.neighbors[i][j];
      vector<unsigned> const & nb = fg.neighbors[other];
      BOOST_CHECK( find( nb.begin(), nb.end(), i) < nb.end() ); 
    }
  
  // check factor indexing and neighbors based on known links
  vector<unsigned> nb = fg.neighbors[ fg.factors[0] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.variables[0] ) < nb.end() ); // fac 0 links var 0

  nb = fg.neighbors[ fg.factors[2] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.variables[0] ) < nb.end() ); // fac 2 links var 0

  nb = fg.neighbors[ fg.factors[1] ];
  BOOST_CHECK( not ( find( nb.begin(), nb.end(), fg.variables[1] ) < nb.end() ) ); // fac 1 does not link var 1

  // check variable indexing and neighbors based on known links
  nb = fg.neighbors[ fg.variables[0] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.factors[0] ) < nb.end() ); // var 0 links fac 0

  nb = fg.neighbors[ fg.variables[2] ];
  BOOST_CHECK( find( nb.begin(), nb.end(), fg.factors[1] ) < nb.end() ); // var 2 links fac 1

  nb = fg.neighbors[ fg.variables[1] ];
  BOOST_CHECK( not ( find( nb.begin(), nb.end(), fg.factors[1] ) < nb.end() ) ); // var 1 does not link fac 1

  //check potentials
  BOOST_CHECK( fg.getFactor(0).getPotential().size1() == 4 );
  BOOST_CHECK( fg.getFactor(0).getPotential().size2() == 4 );
  BOOST_CHECK( fg.getFactor(1).getPotential().size1() == 4 );
  BOOST_CHECK( fg.getFactor(1).getPotential().size2() == 4 );
  BOOST_CHECK( fg.getFactor(2).getPotential().size1() == 1 );
  BOOST_CHECK( fg.getFactor(2).getPotential().size2() == 4 );
}

BOOST_AUTO_TEST_CASE(initMessages_1) 
{
  // TESTS IMPLEMENTATION
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  fg.initMessages();

  unsigned i = 0;
  unsigned j = 1;
  unsigned nb = fg.neighbors[i][j];
  unsigned k  = getIndex(fg.neighbors[nb], i);
  BOOST_CHECK(fg.inMessages_[i][j]  == & fg.outMessages_[nb][k]);
  BOOST_CHECK(fg.inMessages_[nb][k] == & fg.outMessages_[i][j]);
}

BOOST_AUTO_TEST_CASE(sumProduct_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  stateMaskVec_t stateMasks(3);

  // Test 1
  number_t p = fg.calcNormConst(stateMasks);
  BOOST_CHECK_CLOSE(p, 1.0, EPS);

  // Test 2
  stateMask_t v0 = stateMask_t(4, 0);
  v0(0) = 1;
  stateMasks.at(0) = &v0;
  p = fg.calcNormConst(stateMasks);
  BOOST_CHECK_CLOSE( p, 0.25, EPS);
}

BOOST_AUTO_TEST_CASE(calcVariableMarginals_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  stateMaskVec_t stateMasks(3);
  stateMask_t v2 = stateMask_t(4, 0);
  v2(0) = 1;
  stateMasks.at(2) = &v2;
  fg.runSumProduct(stateMasks);

  // calcVariableMarginals
  vector<vector_t> variableMarginals;
  fg.calcVariableMarginals(variableMarginals, stateMasks);

  // test that marginals sum to one
  for (unsigned i = 0; i < variableMarginals.size(); i ++)
    BOOST_CHECK_CLOSE( sum(variableMarginals[i]), 1.0, EPS );

  // test specific marginals
  BOOST_CHECK_CLOSE( variableMarginals[0](0), 0.7, EPS);
}

BOOST_AUTO_TEST_CASE(maxSum_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  
  stateMaskVec_t stateMasks(3);
  stateMask_t v0 = stateMask_t(4, 0), v1 = stateMask_t(4,0), v2 = stateMask_t(4,0);
  v0(0) = 1;
  v1(0) = 1;
  v2(1) = 1;
  stateMasks.at(0) = &v0;
  stateMasks.at(1) = &v1;
  stateMasks.at(2) = &v2;

  // run max sum
  vector<unsigned> maxVar;
  fg.initMaxVariables(maxVar);
  number_t p = fg.runMaxSum(stateMasks, maxVar);

  BOOST_CHECK_CLOSE(p, 0.25*0.7*0.1, EPS);
}

BOOST_AUTO_TEST_CASE(calcFactorMarginals_1) 
{
  // setup
  DFG fg = mkSimpleFactorGraph(); // fg defined a two leaf phylo tree (see function definition)
  stateMaskVec_t stateMasks(3);
  stateMask_t v2 = stateMask_t(4, 0);
  v2(0) = 1;
  stateMasks.at(2) = &v2;

  fg.runSumProduct(stateMasks);

  // calcFactorMarginals
  vector<matrix_t> factorMarginals;
  fg.calcFactorMarginals(factorMarginals);

  // test that marginals sum to one
  for (unsigned i = 0; i < factorMarginals.size(); i ++)
    BOOST_CHECK_CLOSE( sumMatrix(factorMarginals[i]), 1.0, EPS );

  // test specific entries
  BOOST_CHECK_CLOSE( factorMarginals[1](0,0), 0.7, EPS);
  BOOST_CHECK_CLOSE( factorMarginals[1](1,0), 0.1, EPS);
  BOOST_CHECK_CLOSE( factorMarginals[1](2,0), 0.1, EPS);
  BOOST_CHECK_CLOSE( factorMarginals[1](3,0), 0.1, EPS);
}
