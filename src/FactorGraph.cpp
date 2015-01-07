/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "FactorGraph.h"

// debug
#include <boost/numeric/ublas/io.hpp>

namespace phy {

  FGBaseNode::FGBaseNode(vector<FGVertex_t> const & neighbors, vector<unsigned> const & dimension, bool isFactor) : neighbors(neighbors), inMessages(neighbors.size(), NULL), dimension(dimension), received( neighbors.size(), false), isFactor(isFactor)
  {
    if (isFactor) 
      BOOST_FOREACH(unsigned dim, dimension)
	outMessages.push_back( vector_t(dim) );
    else 
      for (unsigned i = 0; i < neighbors.size(); i++)
	outMessages.push_back( vector_t(dimension[0]) );

    // init to zero
    BOOST_FOREACH(vector_t & v, outMessages) 
      v.clear();
  }

  void DiscreteVariableNode::calcOutSumMessage(unsigned const nbIndex)
  {
    unsigned n = neighbors.size();
    assert (n > 0);

    vector_t & v = outMessages[ nbIndex ];

    if (n == 1) {
      for (unsigned i = 0; i < v.size(); i++)
	v[i] = 1;
    }

    else if (n == 2) {
      for (unsigned i = 0; i < n; i++)
	if (i != nbIndex)
	  v = *inMessages[i];
    }

    else { // n > 2
      unsigned i = 0;
      for (; i < n; i++)
	if (i != nbIndex)
	  break;
      v = *inMessages[i];
      i++;
      for (; i < n; i++)
	if (i != nbIndex) 
	  for (unsigned j = 0; j < v.size(); j++)
	    v[j] *= (*inMessages[i])[j];
    }

    // observed variable?
    if (observed)
      for (unsigned i = 0; i < v.size(); i++) {
	v[i] *= (*stateMask)[i];
      }
  }


  void DiscreteFactor2DNode::calcOutSumMessage(unsigned nbIndex) 
  {
    if (nbIndex == 0)
      outMessages[nbIndex] = prod(*inMessages[1], potential);
    else if(nbIndex == 1)
      outMessages[nbIndex] = prod(potential, *inMessages[0]);
    else {
      errorAbort("BaseFactor2DNode::calcOutSumMessage: nbIndex out of range.");
    }
  }


  void inwardCalcSumMessage(vector<FGBaseNode *> nodes, unsigned root)
  {
    inwardCalcSumMessageRec(nodes, root, root);
  }

  // initiate recursion with: inwardsSumMessageRec(nodes, root, root)
  vector_t const * inwardCalcSumMessageRec(vector<FGBaseNode *> nodes, unsigned ndIdx, unsigned sender)
  {
    FGBaseNode & nd(* nodes[ndIdx]);

    // recursively fill out inMessages
    unsigned n = nd.neighbors.size();
    for (unsigned i = 0; i < n; i++)
      if (nd.neighbors[i] != sender) 
	nd.inMessages[i] = inwardCalcSumMessageRec(nodes, nd.neighbors[i], ndIdx);

    // calc and return outMessage
    if (ndIdx == sender) // this is root node of recursion
      return NULL;
    else {
      unsigned nbId = nd.getNeighborIndex(sender);
      nd.calcOutSumMessage(nbId);
      return & nd.outMessages[nbId]; 
    }
  }
 
  void outwardCalcSumMessage(vector<FGBaseNode *> nodes, unsigned root)
  {
    outwardCalcSumMessageRec(nodes, root, root, NULL);
  }

  // inwardCalcSumMessageRec must be run first
  void outwardCalcSumMessageRec(vector<FGBaseNode *> nodes, unsigned ndIdx, unsigned sender, vector_t const * message)
  {
    FGBaseNode & nd(* nodes[ndIdx]);
    if (ndIdx != sender)
      nd.inMessages[ nd.getNeighborIndex(sender) ] = message;

    unsigned n = nd.neighbors.size();
    for (unsigned i = 0; i < n; i++)
      if (nd.neighbors[i] != sender) {
	nd.calcOutSumMessage(i);
	outwardCalcSumMessageRec(nodes, nd.neighbors[i], ndIdx, & nd.outMessages[i]);
      }
  }


  // fill out all inMessages of all nodes
  void calcAllSumMessages(vector<FGBaseNode *> nodes, unsigned root)
  {
    inwardCalcSumMessageRec(nodes, root, root);
    outwardCalcSumMessageRec(nodes, root, root, NULL);
  }

} // end namespace phy
