/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __FactorGraph_h
#define __FactorGraph_h

#include "Tree.h"
#include "PhyDef.h"
#include "utils.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;

// bool initPotentials
// bool initObservedVariables;

namespace phy {

typedef ublas::vector<number_t> vector_t;
typedef ublas::matrix<number_t> matrix_t;

  /** Graph type for factor grahs. Note that factor graphs must be trees. */
  typedef adjacency_list <vecS, vecS, undirectedS> FG_t;

  /** Traits type for convenience */
  typedef graph_traits<FG_t> FGTraits_t;

  /** Vertex type */
  typedef FG_t::vertex_descriptor FGVertex_t;

  /** Edge type */
  typedef FG_t::edge_descriptor FGEdge_t;

  
  struct FGBaseNode {
    virtual ~FGBaseNode() {};

    FGBaseNode(vector<FGVertex_t> const & neighbors, vector<unsigned> const & dimension, bool isFactor);

    /** Calc outgoing sum message for neighbor with index nbIndex and store it in outMessages */
    virtual void calcOutSumMessage(unsigned nbIndex) = 0;

    unsigned getNeighborIndex(unsigned id) 
    {
      unsigned i = 0;
      for (; i < neighbors.size(); i++)
	if (neighbors[i] == id)
	  return i;
      errorAbort("From getNeighborIndex: requested index not found.");
      return 0; // will never reach this
    }

    /** Vector of neighbors. */
    vector<FGVertex_t> const neighbors;

    /** vector of outgoing messages */
    vector<vector_t> outMessages;

    /** vector of incoming messages */
    vector<vector_t const *> inMessages;

    /** Vector of dimensions. Zero, by definition, denotes infinite
	dimensionality and should be used for continuous
	variables. For discrete variable nodes the sole value defines
	the dimensionality of the random variable. For factor nodes
	the vector defines the dimensionality of the potential. */ 
    vector<unsigned> const dimension;

    /** received[i] is true if message from neighbors[i] has been received. */
    vector<bool> received;

    /** True if node is a factor node and false if it is a variable node */
    bool isFactor;
  };

  struct DiscreteVariableNode : public FGBaseNode
  {
    DiscreteVariableNode(vector<FGVertex_t> const & neighbors, vector<unsigned> const & dimension) : FGBaseNode(neighbors, dimension, false), observed(false) {}

    void calcOutSumMessage(unsigned nbIndex);

    /** Observed variable */
    bool observed;

    /** Observation mask */
    ublas::vector<bool> const * stateMask;
  };
 
  struct DiscreteFactor1DNode : public FGBaseNode
  {
  public:
    DiscreteFactor1DNode(vector<FGVertex_t> const & neighbors, 
		     vector<unsigned> const & dimension, 
		     vector_t const & potential) 
      : FGBaseNode(neighbors, dimension, true), potential(potential) {}

    void calcOutSumMessage(unsigned nbIndex) {outMessages[nbIndex] = potential;}

  protected:

    vector_t potential;
  };

  struct DiscreteFactor2DNode : public FGBaseNode
  {
  public:
    DiscreteFactor2DNode(vector<FGVertex_t> const & neighbors, 
			 vector<unsigned> const & dimension, 
			 matrix_t const & pot) 
      : FGBaseNode(neighbors, dimension, true), potential(pot) {}

    void calcOutSumMessage(unsigned nbIndex);

  protected:

    /** Defines the two dimensional potential f(x_i, x_j). The order of the variables must be the same as the order given in the neighbors field (defined in FGBaseNode). */
    matrix_t potential;
  };


  // free functions

  /** Initiates recursive call to the sum-product algorithm and fill out all in-messages from the leaves to the root.*/
  void inwardCalcSumMessage(vector<FGBaseNode *> nodes, unsigned root);

  /** Calculate and send the messages of the sum-product algorithm recursively from the leaves to the root. Messages are stores in the nodes of the nodes vector.*/
  vector_t const * inwardCalcSumMessageRec(vector<FGBaseNode *> nodes, unsigned ndIdx, unsigned sender);

  /** Initiates recursive call to the sum-product algorithm and out all in-messages from the root to the leaves. Precondition: call inwardCalcSumMessage. */
  void outwardCalcSumMessage(vector<FGBaseNode *> nodes, unsigned root);

  /** Calculate and send the messages of the sum-product algorithm recursively from the root to the leaves (inwardCalcSumMessageRec must first be run). Messages are stores in the nodes of the nodes vector. */
  void outwardCalcSumMessageRec(vector<FGBaseNode *> nodes, unsigned ndIdx, unsigned sender, vector_t const * message);

  /** Run the sum-product algorithm with message passing in both directions along all edges of the graph. Messages are stores in the nodes of the nodes vector. */
  void calcAllSumMessages(vector<FGBaseNode *> nodes, unsigned root);

} // end namespace phy





#endif  //__FactorGraph_h
