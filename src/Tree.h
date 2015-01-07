/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Tree_h
#define __Tree_h

#include <iostream>
#include <utility>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_traits.hpp>
#include "PhyDef.h"

using namespace std;
using namespace boost;

namespace phy {

  /** Structure that holds information of a phylogenetic tree. */
  struct TreeInfo {
    unsigned root;               // node index of root
    vector<string> nodeMap;      // defines names of nodes, which may be the empty string, and defines an enumeration (indexing) of nodes.
    vector<number_t> branchMap;  // defines branch lengths and branch enumeration.
    vector<vector<unsigned> > branchNeighbors; // defines node neighbors of branches. I.e., branchNeighbors[i] gives the indexes of the two nodes connected by branch i. The two neighbors of a branch are ordered by proximity to the root, i.e., parent node followed by child node.
  };


  /** Tree classes, based on boost graph library, and accompanying
   * tree traversal functions follows below. These are no longer used
   * in the factor graph framework. */

  /** Bundled internal Vertex properties*/
  struct VertexProp { 
    string name;
  };
	
  /** Bundled internal Edge properties*/
  struct EdgeProp { 
    EdgeProp(double l = -1.0, int i = -1) : length(l), id(i) {};
    double length;
    /** id can be used as index for a vector of external data */ 
    int id;
  };

  /** Graph type for trees */
  typedef adjacency_list <vecS, vecS, bidirectionalS, VertexProp, EdgeProp, boost::no_property, vecS> tree_t;

  /** Traits type for convenience */
  typedef graph_traits<tree_t> treeTraits_t;

  /** Vertex type */
  typedef tree_t::vertex_descriptor vertex_t;

  /** Edge type */
  typedef tree_t::edge_descriptor edge_t;

  /** Generic function for doing preorder tree traversals. preVisit is
      a functor, which is called on every node before visiting any
      child nodes. */
  template <class PreVisitor>
  void treeTraversal(tree_t const & tree, vertex_t const & v, PreVisitor & preVisit)
  {
    preVisit(v);
    treeTraits_t::adjacency_iterator ai, ai_end;
    for ( tie(ai, ai_end) = adjacent_vertices(v, tree); ai != ai_end; ++ai)
      treeTraversal(tree, *ai, preVisit);
  }

  /** Non-const-tree version of the above. */
  template <class PreVisitor>
  void treeTraversal(tree_t & tree, vertex_t const & v, PreVisitor & preVisit)
  {
    preVisit(v);
    treeTraits_t::adjacency_iterator ai, ai_end;
    for ( tie(ai, ai_end) = adjacent_vertices(v, tree); ai != ai_end; ++ai)
      treeTraversal(tree, *ai, preVisit);
  }

  /** Generic function for doing pre-order tree and post-order
      traversals. preVisit and postVisit are functors, which are
      called bafore and after child vertices have been visited,
      respectively. */
  template <class PreVisitor, class PostVisitor>
  void treeTraversal(tree_t const & tree, vertex_t const & v, PreVisitor & preVisit, PostVisitor & postVisit)
  {
    preVisit(v);
    treeTraits_t::adjacency_iterator ai, ai_end;
    for ( tie(ai, ai_end) = adjacent_vertices(v, tree); ai != ai_end; ++ai)
      treeTraversal(tree, *ai, preVisit, postVisit);
    postVisit(v);
  }

  /** Same as the above with the addition of envoking inVisit on nodes in-between visiting any children */
  /** note that inVisitor visits nodes in-between visiting children. It visits a node n times, where n equals the number of children minus one */
  template <class PreVisitor, class PostVisitor, class InVisitor>
  void treeTraversal(tree_t const & tree, vertex_t const &v, PreVisitor & preVisit, PostVisitor & postVisit, InVisitor & inVisit)
  {
    preVisit(v);
    treeTraits_t::adjacency_iterator ai, ai_end;
    for ( tie(ai, ai_end) = adjacent_vertices(v, tree); ai != ai_end; ++ai) {
      treeTraversal(tree, *ai, preVisit, postVisit, inVisit);
      if (ai + 1 != ai_end)
	inVisit(v);
    }
    postVisit(v);
  }

} // end namespace phy


#endif  //__Tree_h
