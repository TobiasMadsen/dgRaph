/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __DiscreteFactorGraph_h
#define __DiscreteFactorGraph_h

#include "PhyDef.h"
#include "utils.h"
#include "utilsLinAlg.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <utility>
#include <algorithm>

namespace phy {

using namespace std;

  struct DFGNode
  {
    /**Variable node constructor*/
    DFGNode(unsigned dimension);
    
    /**Factor node constructor*/
    DFGNode(xmatrix_t const & potential);

    /**Factor node with exp*/
    DFGNode(xmatrix_t const & potential, xmatrix_t const & fun_a_init, xmatrix_t const & fun_b_init);

    bool isFactor;      // true if factor node, false if variable node
    unsigned dimension; // dimension of variable or dimension of potential
    xmatrix_t potential; // n * m dimensional matrix defining factor potential. If n == 1, then the potential is treated as one-dimensional. n = m = 0 for variable nodes.
    xmatrix_t fun_a; //See sumproduct note
    xmatrix_t fun_b; // Expectancies 
  };

  class DFG 
  {
  public: 

    DFG(vector<unsigned> const & varDimensions, // dimensions of random variables. Their enumeration is implicit
	vector<xmatrix_t> const & facPotentials, // factor potentials (see definition in DFGNode)
	vector<vector<unsigned> > const & facNeighbors);  // see below
    // facNeighbors[i] (nbs) defines the node neighbors of the i'th
    // factor. Currently, factors can have at most two
    // neighbors. (nbs[0] should refer to the parent node in the case
    // of directed graphical models.) In general the dimension of the
    // nbs[0] node should match the first dimension of the factor
    // potential and the dimension of nbs[1] should match the second
    // dimension of the factor potential. However, if factors only
    // have one neighbor, then the dimension of the nbs[0] node should
    // match the second dimention of the fator potential (i.e. the
    // potential matrix is a 1xn matrix, where n is the dimension of
    // the nbs[0] node). This is checked by consistencyCheck(), which is
    // called by the constructor.
    DFG(vector<unsigned> const & varDimensions,
	vector<xmatrix_t> const & facPotentials,
	vector<vector<unsigned> > const & facNeighbors,
	vector<xmatrix_t> const & facFunA,
	vector<xmatrix_t> const & facFunB);
    

    // same as above but allows factor potentials to be defined as matrix_t types.
    DFG(vector<unsigned> const & varDimensions, // dimensions of random variables. Their enumeration is implicit
	vector<matrix_t> const & facPotentials, // factor potentials (see definition in DFGNode)
	vector<vector<unsigned> > const & facNeighbors);  // see below

    DFG(vector<unsigned> const & varDimensions,
	vector<matrix_t> const & facPotentials,
	vector<vector<unsigned> > const & facNeighbors,
	vector<matrix_t> const & facFunA,
	vector<matrix_t> const & facFunB);

    /** The copy constructor is defined explicitly to ensure the
	inMessages_ pointers are wiped out on copy. This means that
	the states of the dynammic programming tables are not kept on
	copy.*/
    DFG(DFG const & rhs) : nodes(rhs.nodes), neighbors(rhs.neighbors), variables(rhs.variables), factors(rhs.factors), components(rhs.components), roots(rhs.roots) {}

    // public data
    vector<DFGNode> nodes;       // Contains and enumerates all the nodes of the graph.
    vector< vector< unsigned > > neighbors; // Defines graph structure in terms of neighboring nodes 
    vector<unsigned> variables;  // enumerates variable nodes
    vector<unsigned> factors;    // enumerates factor nodes
    vector<unsigned> components;  // assigns each node to a connected component 0 indicates no assignment
    vector<unsigned> roots; // a root node in each component

    // Read-only accessors
    const vector<xmatrix_t> & getXFactorMarginals(){ return xFactorMarginals;}
    const vector<xvector_t> & getXVariableMarginals(){ return xVariableMarginals;}
    const vector<matrix_t> & getFactorMarginals(){ return mFactorMarginals;}
    const vector<vector_t> & getVariableMarginals(){ return mVariableMarginals;}

    /** Reset specific or all factor potentials. Useful with factor optimization, e.g., by EM. */
    // xmatrix_t versions
    void resetFactorPotential(xmatrix_t const & facPotVec, unsigned facId);
    void resetFactorPotentials(vector<xmatrix_t> const & facPotVecSubSet, vector<unsigned> const & facMap);  // facMap maps positions of facPotVecSubSet onto facPotVec
    void resetFactorPotentials(vector<xmatrix_t> const & facPotVec);
    
    // overloaded matrix_t versions
    void resetFactorPotential(matrix_t const & facPotVec, unsigned facId);
    void resetFactorPotentials(vector<matrix_t> const & facPotVecSubSet, vector<unsigned> const & facMap);  // facMap maps positions of facPotVecSubSet onto facPotVec
    void resetFactorPotentials(vector<matrix_t> const & facPotVec);

    /** Calc the normalization constant (Z) for the factor graph: Z =
	sum p(x_root). Runs a single pass of the sum-product algorithm
	(half the computations of the full sumProduct algorithm needed
	for calculating all local marginals). Not sufficient if
	variable marginals are also needed.  */
    xnumber_t calcNormConst(stateMaskVec_t const & stateMasks);
    xnumber_t calcNormConst(stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
 
    /** Run the sum product algorithm. Postcondition: all inMessages
	and outMessages are computed. The stateMask: Observed data is
	given through the stateMask. If no data is observed for
	variable i, then stateMasks[i] is a NULL pointer, otherwise
	stateMasks[i] points to a stateMask with value one if a state
	is observed and zero otherwise. */
    void runSumProduct(stateMaskVec_t const & stateMasks);  
    void runSumProduct(stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;

    /** Precondition: runSumProduct has been called (setting all out/inMessages). Calc the normalization constant (Z). */
    xnumber_t calcNormConst2(stateMaskVec_t const & stateMasks) const;
    xnumber_t calcNormConst2(stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages) const;

    /** Precondition: runSumProcut has been called (setting all out/inMessages). If variableMarginals is an empty vector ::initVariableMarginals will be called*/
    void calcVariableMarginals(stateMaskVec_t const & stateMasks);
    void calcVariableMarginals(vector<xvector_t> & variableMarginals, stateMaskVec_t const & stateMasks);
    void calcVariableMarginals(vector<xvector_t> & variableMarginals, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages) const;

    /** Precondition: runSumProcut has been called (setting all out/inMessages). If factorMarginals is an empty vector ::initFactorMarginals will be called*/
    void calcFactorMarginals();
    void calcFactorMarginals(vector<xmatrix_t> & factorMarginals);
    void calcFactorMarginals(vector<xmatrix_t> & factorMarginals, vector<vector<xvector_t const *> > & inMessages) const;

    /** Run the max sum algorithm (should be called max product, since
	computation is not done in log space as is usually the
	case). Also performs the backtrack stage of the
	algorithm. Returns the probability of the most probable
	outcome. Postcondition: maxVariable defines the most probable
	outcome, maxVariable[i] is the state of variable i. */
    xnumber_t runMaxSum(stateMaskVec_t const & stateMasks, vector<unsigned> & maxVariables);
    xnumber_t runMaxSum(stateMaskVec_t const & stateMasks, vector<unsigned> & maxVariables, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const; 

    /** Calculate likelihood for a full observation */
    number_t calcFullLikelihood( vector<unsigned> const & sample);

    /** Precondition: runSumProduct has been called */
    xnumber_t calcExpect(); ///< Using naive algorithm
    xnumber_t calcExpect2(stateMaskVec_t const & stateMasks); ///< Using (if implemented correctly) memory efficient algorithm

    /** No precondition give the two kinds of potentials*/
    pair<xnumber_t,xnumber_t> calcExpectancies(vector<xmatrix_t> const & fun_a, vector<xmatrix_t> const & fun_b, stateMaskVec_t const & stateMasks);

    /** write factor graph in dot format (convert to ps using: cat out.dot | dot -Tps -o out.ps ) */
    void writeDot(string const & fileName);

    // convenience functions
    DFGNode const & getFactor(unsigned facId) const {return nodes[ factors[ facId ] ];} 
    DFGNode const & getVariable(unsigned varId) const {return nodes[ variables[ varId ] ];}
    unsigned convNodeToFac(unsigned ndId) const {return getIndex(factors, ndId);} 
    unsigned convNodeToVar(unsigned ndId) const {return ndId;} // exploiting how node ids were constructed [getIndex(variables, ndId);}] 
    unsigned convFacToNode(unsigned facId) const {return factors[facId];} 
    unsigned convVarToNode(unsigned varId) const {return variables[varId];} 
    vector<unsigned> const & getFactorNeighbors(unsigned facId) const {return neighbors[ factors[ facId ] ];}     // neighbors are defined in terms of 'nodes' indices
    vector<unsigned> const & getVariableNeighbors(unsigned varId) const {return neighbors[ variables[ varId ] ];} // neighbors are defined in terms of 'nodes' indices

    // functions aiding in setting up data structures. Useful if external data structures are used.
    /** Initialize outMessages to be of the right size and inMessages to point to the corresponding out messages */
    void initMessages(vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
    void initVariableMarginals(vector<xvector_t> & variableMarginals) const;
    void initFactorMarginals(vector<xmatrix_t> & factorMarginals) const;
    /** Sets up the data structure needed for maxSum. Indexing:
	maxNeighborStates[i][j] returns a vector of unsigned (states)
	corresponding to the j'th neigbor of the i'th factor. */
    void initMaxNeighbourStates(vector<vector<vector<unsigned> > > & maxNeighborStates) const;
    void initMaxVariables(vector<unsigned> & maxVariables) const;


    /** Consistency check of factor dimensions versus number of neighbors and of potential dimensions and neighbors (var) dimensions */
    void consistencyCheck();

    /** Make a sample from conditional distribution induced by factorgraph, precondition sumProduct has been run */
    void sample(boost::mt19937 & gen, vector<vector_t> const & varMarginals, vector<matrix_t> const & facMarginals, vector<unsigned> & sample);
    vector<unsigned> sample(boost::mt19937 & gen, vector<vector_t> const & varMarginals, vector<matrix_t> const & facMarginals);
    void sampleIS(boost::mt19937 & gen, vector<vector_t> const & varMarginals, vector<matrix_t> const & facMarginals, vector<vector_t> const & ISVarMarginals, vector<matrix_t> const & ISFacMarginals, vector<unsigned> & sim, number_t & weight);


    /** Precondition: calcFactorMarginals and calcVariableMarginals has been called on members factorMarginals and variableMarginals */
    void simulateVariable(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<unsigned> & sim);
    void simulateFactor(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<unsigned> & sim);
    void simulateVariableIS(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<matrix_t> const & ISFacMarginals, vector<unsigned> & sim, number_t & weight);
    void simulateFactorIS(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<matrix_t> const & ISFacMarginals,  vector<unsigned> & sim, number_t & weight);
    
    /** write factor info to str */
    void writeInfo( ostream & str, vector<string> const & varNames = vector<string>(), vector<string> const & facNames = vector<string>() );

  protected:

    // initialization called by constructors
    void init(vector<unsigned> const & varDimensions, vector<xmatrix_t> const & facPotentials, vector<vector<unsigned> > const & facNeighbors);

    void init(vector<unsigned> const & varDimensions, vector<xmatrix_t> const & facPotentials, vector<vector<unsigned> > const & facNeighbors, vector<xmatrix_t> const & facFunA, vector<xmatrix_t> const & facFunB);

    // return string with info on factor i
    string factorInfoStr( unsigned const i, vector<string> varNames = vector<string>(), vector<string> facNames = vector<string>() );
    string variableInfoStr( unsigned const i, vector<string> varNames = vector<string>(), vector<string> facNames = vector<string>() );

    // helper functions for sumProduct();
    void runSumProductInwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
    void runSumProductOutwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
    void calcSumProductMessageFactor(unsigned current, unsigned receiver, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
    void calcSumProductMessageFactor(unsigned current, unsigned receiver, vector<xvector_t const *> const & inMes, xvector_t & outMes) const;
    void calcSumProductMessageVariable(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
    void calcSumProductMessageVariable(unsigned current, unsigned receiver, stateMask_t const * stateMask, vector<xvector_t const *> const & inMes, xvector_t & outMes) const;
    void calcSumProductMessage(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const;
    xnumber_t calcNormConst(unsigned varId, stateMask_t const * stateMask, vector<xvector_t const *> const & inMes) const;

    // helper functions for maxSum();
    unsigned maxNeighborDimension(vector<unsigned> const & nbs) const;
    void runMaxSumInwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const;
    void calcMaxSumMessage(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const;
    void calcMaxSumMessageFactor(unsigned current, unsigned receiver, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const;
    void calcMaxSumMessageFactor(unsigned current, unsigned receiver, vector<xvector_t const *> const & inMes, xvector_t & outMes, vector<vector<unsigned> > & maxNBStates) const;
    void backtrackMaxSumOutwardsRec(vector<unsigned> & maxVariables, unsigned current, unsigned sender, unsigned maxState, vector<vector<vector<unsigned> > > & maxNeighborStates) const;

    //helper functions for calcExpect2
    void runExpectInwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages2, vector<vector<xvector_t> > & outMessages2) const;
    void runExpectOutwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages2, vector<vector<xvector_t> > & outMessages2) const;
    void calcExpectMessageFactor(unsigned current, unsigned receiver, vector<vector<xvector_t const *> > & inMessages2, vector<vector<xvector_t> > & outMessages2) const;
    void calcExpectMessageFactor(unsigned current, unsigned receiver, vector<xvector_t const *> const & inMes, xvector_t & outMes2, vector<xvector_t const *> const & inMes2) const;
    void calcExpectMessageVariable(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages2, vector<vector<xvector_t> > & outMessages2) const;
    void calcExpectMessageVariable(unsigned current, unsigned receiver, stateMask_t const * stateMask, vector<xvector_t const *> const & inMes2, xvector_t & outMes2, vector<xvector_t const *> const & inMes) const;
    void calcExpectMessage(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages2, vector<vector<xvector_t> > & outMessages2) const;


    //Functions for calculating expectancies
    //See note: sumProduct.pdf
    void runExpectanciesInwardsRec(unsigned current, unsigned sender, vector<xmatrix_t> const & fun_a, vector<xmatrix_t> const & fun_b, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMu, vector<vector<xvector_t> > & outMu, vector<vector<xvector_t const *> > & inLambda, vector<vector<xvector_t> > & outLambda) const;
//void runExpectanciesOutwardsRec(unsigned current, unsigned sender, vector<xmatrix_t> const & fun_a, vector<xmatrix_t> const & fun_b, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMu, vector<vector<xvector_t> > & outMu, vector<vector<xvector_t const *> > & inLambda, vector<vector<xvector_t> > & outLambda) const;
    void calcExpectanciesMessageFactor(unsigned current, unsigned receiver, vector<xmatrix_t> const & fun_a, vector<xmatrix_t> const & fun_b, vector<vector<xvector_t const *> > & inMu, vector<vector<xvector_t> > & outMu, vector<vector<xvector_t const *> > & inLambda, vector<vector<xvector_t> > & outLambda) const;
    void calcExpectanciesMessageFactor(unsigned current, unsigned receiver, vector<xmatrix_t> const & fun_a, vector<xmatrix_t> const & fun_b, vector<xvector_t const *> const & inMesMu, xvector_t & outMesMu, vector<xvector_t const *> const & inMesLambda, xvector_t & outMesLambda) const;
    void calcExpectanciesMessageVariable(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMu, vector<vector<xvector_t> > & outMu, vector<vector<xvector_t const *> > & inLambda, vector<vector<xvector_t> > & outLambda) const;
    void calcExpectanciesMessageVariable(unsigned current, unsigned receiver, stateMask_t const * stateMask, vector<xvector_t const *> const & inMesMu, xvector_t & outMesMu, vector<xvector_t const *> const & inMesLambda, xvector_t & outMesLambda) const;
    void calcExpectanciesMessage(unsigned current, unsigned sender, vector<xmatrix_t> const & fun_a, vector<xmatrix_t> const & fun_b, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMu, vector<vector<xvector_t> > & outMu, vector<vector<xvector_t const *> > & inLambda, vector<vector<xvector_t> > & outLambda) const;


    /** Convert to boost factor graph (only captures graph structure */ 
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> mkBoostGraph();

    /** write factor graph in dot format to string. */
    string writeDot();

    // Initialize data structures, only need to run once
    void initMessages();
    void initMaxNeighbourStates();
    void initComponents();

    // private data
    // convenience data structures -- perhaps make public
    vector<vector<xvector_t const *> > inMessages_;
    vector<vector<xvector_t> > outMessages_;
    vector<vector<vector<unsigned> > > maxNeighborStates_;

    vector<vector<xvector_t const *> > inMessages2_; //For messages of the second type
    vector<vector<xvector_t> > outMessages2_;

    //Structures for calculation of expectancies
    //See note sumproduct.pdf
    vector<vector<xvector_t const *> > inMu_;
    vector<vector<xvector_t> > outMu_;
    vector<vector<xvector_t const *> > inLambda_;
    vector<vector<xvector_t> > outLambda_;

    vector<vector_t> mVariableMarginals;
    vector<matrix_t> mFactorMarginals;
    vector<xvector_t> xVariableMarginals;
    vector<xmatrix_t> xFactorMarginals;
  };

  // free functions for dealing vith a vector of observation sets

  /** Calculates the normalization constant for a vector of stateMasks
      (observations). If the factors represent proper probability
      distributions, the norm constant equal the likelihood of the
      observations under the model */
  vector<xnumber_t> calcNormConsMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg);
  void calcNormConsMultObs(vector<xnumber_t> & result, stateMask2DVec_t const & stateMask2DVec, DFG & dfg);

  /** Calculates the accumulated variable marginals over all
      stateMask vectors. Useful for EM based parameter estimation,
      etc.. */
  vector<vector_t> calcVarAccMarMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg);
  /** Result is assumed of correct size (use dfg.initVariableMarginals above).  The accumulated counts will be added, therefore reset to zero as needed (use resetVarMarVec below).*/
  void calcVarAccMarMultObs(vector<vector_t> & result, stateMask2DVec_t const & stateMask2DVec, DFG & dfg);

  /** Calculates the accumulated factor marginals over all
      stateMask vectors. Useful for EM based parameter estimation,
      etc.. */
  vector<matrix_t> calcFacAccMarMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg);
  /** Result is assumed of correct size (use dfg.initFactorMarginals above).  The accumulated counts will be added, therefore reset to zero as needed (use resetFacMarVec below).*/
  void calcFacAccMarMultObs(vector<matrix_t> & result, stateMask2DVec_t const & stateMask2DVec, DFG & dfg);

  /** Calculates the accumulated variable and factor marginals over
      all stateMask vectors. */
  /** varResult and facResult are assumed of correct sizes (use dfg.init* member functions above).  The accumulated counts will be added, therefore reset to zero as needed (use reset*Vec below).*/
  void calcVarAndFacAccMarMultObs(vector<matrix_t> & varResult, vector<matrix_t> & facResult, stateMask2DVec_t const & stateMask2DVec, DFG & dfg);

  /** Calculates the most probable combined state assignment to variables for multiple stateMask (observation) vectors. For ease of use, size of result references will be adjusted is of wrong size. */
  void calcMaxProbStatesMultObs(vector<xnumber_t> & maxProbResult, vector<vector<unsigned> > & maxVarResult, stateMask2DVec_t const & stateMask2DVec, DFG & dfg);
  pair<vector<xnumber_t>, vector<vector<unsigned> > > calcMaxProbStatesMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg);

  /** Calculates the accumulated factor marginals as well as the
      normalization constants (Z) for multiple (second version) or
      single stateMask vectors. The accumulated counts will be added,
      therefore reset to zero as needed (use reset*Vec below). Useful
      for EM based parameter estimation. */
  // tmpFacMar is used as workspace. accFacMar must be of correct size (use dfg.initFactorMarginals(tmpFacMar) )
  void calcFacAccMarAndNormConst(vector<matrix_t> & accFacMar, vector<xmatrix_t> & tmpFacMar, xnumber_t & normConst, stateMaskVec_t const & stateMaskVec, DFG & dfg);
  void calcFacAccMarAndNormConstMultObs(vector<matrix_t> & accFacMar, xvector_t & normConstVec, stateMask2DVec_t const & stateMask2DVec, DFG & dfg);

  /* Initiate data structure for accumulated marginals */
  void initAccVariableMarginals(vector<vector_t> & variableMarginals, DFG const & dfg);
  void initAccFactorMarginals(vector<matrix_t> & factorMarginals, DFG const & dfg);

  ////////////////////////////////////////////////////////////////
  // definition of template and inline functions follows below

  /** helper functions. Interface provided above*/
  template <class T>
  void initGenericVariableMarginals(vector<T> & variableMarginals, DFG const & dfg)
  {
    variableMarginals.clear();
    variableMarginals.resize( dfg.variables.size() );
    for (unsigned i = 0; i < dfg.variables.size(); i++) {
      unsigned dim = dfg.nodes[ dfg.variables[i] ].dimension;
      variableMarginals[i].resize(dim);
    }
    reset(variableMarginals);
  }

  template<class T>
  void initGenericFactorMarginals(vector<T> & factorMarginals, DFG const & dfg)
  {
    factorMarginals.clear();
    factorMarginals.resize( dfg.factors.size() );
    for (unsigned i = 0; i < dfg.factors.size(); i++) {
      unsigned size1 = dfg.nodes[ dfg.factors[i] ].potential.size1();
      unsigned size2 = dfg.nodes[ dfg.factors[i] ].potential.size2();
      factorMarginals[i].resize(size1, size2);
    }
    reset(factorMarginals);
  }

  // the member functions below are made inline for efficiency


  inline void DFG::calcSumProductMessageVariable(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    xvector_t & outMes = outMessages[current][ getIndex( nbs, receiver) ];  // identify message
    vector< xvector_t const *> const & inMes( inMessages[current] );
    stateMask_t const * stateMask = stateMasks[ convNodeToVar(current) ]; 

    calcSumProductMessageVariable(current, receiver, stateMask, inMes, outMes);
  }


  inline void DFG::calcSumProductMessageVariable(unsigned current, unsigned receiver, stateMask_t const * stateMask, vector<xvector_t const *> const & inMes, xvector_t & outMes) const
  {
    vector<unsigned> const & nbs = neighbors[current];

    if (stateMask) {  // observed variable?
      for (unsigned i = 0; i < outMes.size(); i++) 
	outMes[i] = (*stateMask)[i];
    }
    else  // init to 1
      for (unsigned i = 0; i < outMes.size(); i++)
	outMes[i] = 1;
      
    for (unsigned i = 0; i < nbs.size(); i++)
      if (nbs[i] != receiver) 
	for (unsigned j = 0; j < outMes.size(); j++)
	  outMes[j] *= (*inMes[i])[j];
  }


  inline void DFG::calcSumProductMessageFactor(unsigned current, unsigned receiver, vector<vector<xvector_t const *> > & inMessages, vector<vector<xvector_t> > & outMessages) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    xvector_t & outMes = outMessages[current][ getIndex( nbs, receiver) ];  // identify message
    vector< xvector_t const *> const & inMes( inMessages[current] );

    calcSumProductMessageFactor(current, receiver, inMes, outMes);
  }


  inline void DFG::calcSumProductMessageFactor(unsigned current, unsigned receiver, vector<xvector_t const *> const & inMes, xvector_t & outMes) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    DFGNode const & nd = nodes[current];

    // one neighbor 
    if (nd.dimension == 1) {
      for (unsigned i = 0; i < nd.potential.size2(); i++)
	outMes[i] = nd.potential(0, i);
      return;
    }

    // two neighbors  (this is were most time is normally spent in normConst calculations)
    if (nd.dimension == 2) {
      if (nbs[0] == receiver) // factor neighbor closests to root (for directed graphs)
	outMes = prod(nd.potential, *inMes[1]);  // note that *inMes[1] is column-vector (as are all ublas vectors)
      else // nbs[1] == receiver  // factor neighbor furthest away from root (for directed graphs)
	outMes = prod(*inMes[0], nd.potential);
      return;
    }

    // more than two neighbors -- should not happen
    errorAbort("calcSumProductMessageFactor: Error, factor with more than two neighbors");
  }




} // end namespace phy

#endif  //__DiscreteFactorGraph_h
