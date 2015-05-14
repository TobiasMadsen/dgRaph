/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "DiscreteFactorGraph.h"

#include <queue>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace phy {

  DFGNode::DFGNode(unsigned dimension) : 
  isFactor_(false), dimension(dimension) {}


  DFGNode::DFGNode(matrix_t const & potential) : 
    isFactor_(true), potential(potential)
    {
      if (potential.size1() == 1)
	dimension = 1;
      else
	dimension = 2;
    }

  bool DFGNode::isFactor() const{
    return isFactor_;
  }

  unsigned DFGNode::getDimension() const{
    return dimension;
  }

  matrix_t DFGNode::getPotential() const{
    return potential;
  }

  void DFGNode::setPotential(matrix_t const & pot){
    potential = pot;
  }

  // The graph is defined in terms of the factor neighbors. We want
  // links to be represented both ways, i.e., also from var nodes to
  // fac nodes. This function returns the neighbors from the var point
  // of view.

  // precondition: neighbors harbor only one-way links
  void addTwoWayLinks(vector<vector<unsigned> > & neighbors)
  {
    vector<vector<unsigned> > ref(neighbors);
    for (unsigned i = 0; i < ref.size(); i++) 
      for (unsigned j = 0; j < ref[i].size(); j++) {
	unsigned other = ref[i][j];
	neighbors[other].push_back(i);
      }
  }


  DFG::DFG(vector<unsigned> const & varDimensions, 
	   vector<matrix_t> const & facPotentials, 
	   vector<vector<unsigned> > const & facNeighbors)
  {
    init(varDimensions, facPotentials, facNeighbors);
  }

  void DFG::resetFactorPotential(matrix_t const & pot, unsigned facId)
  {
    DFGNode & nd = nodes[ convFacToNode(facId) ];
    if( nd.getPotential().size1() != pot.size1() )
      errorAbort("DFG::resetFactorPotential: nd.getPotential().size1()="+toString(nd.getPotential().size1())+" does not match pot.size1()="+toString(pot.size1()) );
    if( nd.getPotential().size2() != pot.size2() )
      errorAbort("DFG::resetFactorPotential: nd.getPotential().size2()="+toString(nd.getPotential().size2())+" does not match pot.size2()="+toString(pot.size2()) );
    nd.setPotential(pot);
  }


  void DFG::resetFactorPotentials(vector<matrix_t> const & facPotVecSubSet, vector<unsigned> const & facMap)
  {
    assert(facPotVecSubSet.size() == facMap.size() );
    for (unsigned i = 0; i < facMap.size(); i++)
      resetFactorPotential( facPotVecSubSet[i], facMap[i] );
  }


  void DFG::resetFactorPotentials(vector<matrix_t> const & facPotVec)
  {
    assert(facPotVec.size() == factors.size() );
    for (unsigned i = 0; i < facPotVec.size(); i++)
      resetFactorPotential(facPotVec[i], i);
  }

  // Begining of  write dot

  /** functor for use with write_graphviz. Writes out internal properties */
  struct DFGNodeWriter 
  {
    DFGNodeWriter(DFG const & fg) : fg_ (fg) {};

    // writing factor node in dot format
    void mkDotFacNode(std::ostream& out, unsigned const & v)
    {
      out << " [label=\"f" << fg_.convNodeToFac(v) << "\\n" << v << "\", shape=box]" << endl;
    }

    // writing variable node in dot format
    void mkDotVarNode(std::ostream& out, unsigned const & v)
    {
      out << " [label=\"v" << fg_.convNodeToVar(v) << "\\n" << v << "\"]" << endl;
    }
    
    template <class Vertex>
    void operator()(ostream & out, Vertex const & v) 
    {
      if (fg_.nodes[v].isFactor()) 
	mkDotFacNode(out, v);
      else
	mkDotVarNode(out, v);
    }

 
    DFG const & fg_;
  };


  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> DFG::mkBoostGraph()
  {
    set< pair<unsigned, unsigned> > seen;
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> bfg( nodes.size() );
    for (unsigned i = 0; i < neighbors.size(); i++)
      for (unsigned j = 0; j < neighbors[i].size(); j++) {
	unsigned from = i;
	unsigned to = neighbors[i][j];
	if(seen.insert(pair<unsigned, unsigned>(from, to) ).second and seen.insert(pair<unsigned, unsigned>(to, from) ).second ) // link not seen 
	  add_edge(i, neighbors[i][j], bfg);
      }
    return bfg;
  }


  string DFG::writeDot()
  {
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> bfg( mkBoostGraph() );
    stringstream ss;
    write_graphviz(ss, bfg, DFGNodeWriter(*this) );
    return ss.str();
  }


  void DFG::writeDot(string const &fileName)
  {
    ofstream f(fileName.c_str(), ios::out);
    if (!f)
      errorAbort("Cannot open file: " + fileName + "\n");
  
    f << writeDot();
    f.close();
  }

  // end of write dot

  void DFG::writeInfo(ostream & str, vector<string> const & varNames, vector<string> const & facNames)
  {
    str << "Variable info:" << endl << endl;
    for (unsigned i = 0; i < variables.size(); i++)
      str << variableInfoStr( i, varNames, facNames) << endl;
    str << endl;

    str << "Factor info:" << endl << endl;
    for (unsigned i = 0; i < factors.size(); i++)
      str << factorInfoStr(i, varNames, facNames) << endl;
    str << endl;
  }


  string DFG::factorInfoStr(unsigned const i, vector<string> varNames, vector<string> facNames)
  {
    bool useFacNames = (facNames.size() > 0) ? true : false;
    bool useVarNames = (varNames.size() > 0) ? true : false;
    string s = 
      "Factor index: " + toString(i) 
      + ( useFacNames  ? ("\t" + facNames[i] ) : "") + "\n"
      + "Node index: " + toString( convFacToNode(i) ) + "\n"
      + "Variable neighbors: " + toString( getFactorNeighbors(i) )
      + ( useVarNames ? ("\t" + toString( mkSubset( varNames, getFactorNeighbors(i) ) ) ) : "")  + "\n" // exploiting that node and var ids are equal
      + "Dimension: " + toString( getFactor(i).getDimension() ) + "\n"
      + "Potential: ";
    std::stringstream ss;
    ss << getFactor(i).getPotential() << endl;
    s += ss.str();
    return s;
  }


  string DFG::variableInfoStr(unsigned const i, vector<string> varNames, vector<string> facNames)
  {
    bool useFacNames = (facNames.size() > 0) ? true : false;
    bool useVarNames = (varNames.size() > 0) ? true : false;
    vector<unsigned> facSubset;
    if (useFacNames) {
      BOOST_FOREACH(unsigned j, getVariableNeighbors(i) )
	facSubset.push_back( convNodeToFac(j) );
    }
    string s = 
      "Variable index: " + toString(i) 
      + ( useVarNames  ? ("\t" + varNames[i] ) : " ") + "\n"
      + "Node index: " + toString( convVarToNode(i) ) + "\n"
      + "Factor neighbors: " + toString( getVariableNeighbors(i) ) 
      + ( useFacNames ? ("\t" + toString( mkSubset( facNames, facSubset ) ) ) : "")  + "\n" 
      + "Dimension: " + toString( getVariable(i).getDimension() ) + "\n";
    return s;
  }


  void DFG::consistencyCheck()
  {
    // check factor dimension and number of neighbors
    for (unsigned i = 0; i < factors.size(); i++) 
      if (getFactorNeighbors(i).size() != getFactor(i).getDimension())
	errorAbort("Inconsistent graph. \nNumber of neighbors and potential dimenension does not match at factor:\n\n" + factorInfoStr(i) );
    
    // check that factor potential dimensions equal dimension of neigboring var nodes
    for (unsigned i = 0; i < factors.size(); i++) {
      vector<unsigned> nbs = getFactorNeighbors(i);
      unsigned nbCount = nbs.size();
      if (nbCount > 0) {
	unsigned nb = nbs[0];
	if ( ( nbCount != 1 and getFactor(i).getPotential().size1() != getVariable(nb).getDimension() ) or (nbCount == 1 and getFactor(i).getPotential().size2() != getVariable(nb).getDimension()) ) // if ncCount equals 1 then the potential is a prior
	  errorAbort("Inconsistent graph. \nFactor potential dimensions does not match dimension of neighboring variable node (both defined below):\n\n" + factorInfoStr(i)  + "\n" + variableInfoStr(i) + "\n\n" + "getFactor(i).getPotential().size1():\t" + toString(getFactor(i).getPotential().size1()) + "\n");
      }
      if (nbCount > 1) {
	unsigned nb = nbs[1];
	if (getFactor(i).getPotential().size2() != getVariable(nb).getDimension() )
	  errorAbort("Inconsistent graph. \nFactor potential dimensions does not match dimension of neighboring variable node (both defined below):\n\n" + factorInfoStr(i)  + "\n" + variableInfoStr(i) + "getFactor(i).getPotential().size2():\t" + toString(getFactor(i).getPotential().size2()) + "\n" + "getVariable(nb).dimension:\t" + toString(getVariable(nb).getDimension()) + "\n" );
      }
      if (nbCount > 2)
	  errorAbort("Inconsistent graph. \nFactor with more than two neighbors currently not supported:\n" + factorInfoStr(i) );
    }
  }

  void DFG::runSumProduct(stateMaskVec_t const & stateMasks)
  {
    if (inMessages_.size() == 0) 
      initMessages();
    
    runSumProduct(stateMasks, inMessages_, outMessages_);
  }


  void DFG::runSumProduct(stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages) const
  {
    if(roots.size() == 0)
      errorAbort("DiscreFactorGraph.cpp::runSumProduct: No root nodes");
    for(int i = 0; i < roots.size(); ++i){
      unsigned root = roots.at(i);
      runSumProductInwardsRec(root, root, stateMasks, inMessages, outMessages);
      runSumProductOutwardsRec(root, root, stateMasks, inMessages, outMessages);
    }
  }


  void DFG::calcSumProductMessage(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages) const
  {
    if (nodes[current].isFactor())
      calcSumProductMessageFactor(current, receiver, inMessages, outMessages);
    else
      calcSumProductMessageVariable(current, receiver, stateMasks, inMessages, outMessages);
  };


  // send requests from the root out, returns messages from the leaves
  // inwards. Note that messages are passed implicitly through the
  // tables of outMessages (and pointers of inMessages)
  void DFG::runSumProductInwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    // recursively call all nodes
    for (unsigned i = 0; i < nbs.size(); i++) {
      unsigned nb = nbs[i];
      if (nb != sender) 
	runSumProductInwardsRec(nb, current, stateMasks, inMessages, outMessages);
    }

    // calc outMessage and store in outMessages -- inMessages provides links to outMessages, which provides the message passing.
    if (current == sender) // this is root node of recursion
      return;
    calcSumProductMessage(current, sender, stateMasks, inMessages, outMessages); //  receiver = sender
  }


  // send messages from the root outwards to the leaves. Precondition: DFG::runSumProductInwardsRec must have been run on inMessages and outMessages.
  void DFG::runSumProductOutwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    // recursively call all nodes
    for (unsigned i = 0; i < nbs.size(); i++) {
      unsigned nb = nbs[i];
      if (nb != sender) {
	calcSumProductMessage(current, nb, stateMasks, inMessages, outMessages); //  receiver = nb
	runSumProductOutwardsRec(nb, current, stateMasks, inMessages, outMessages);
      }
    }
  }



  number_t DFG::runMaxSum(stateMaskVec_t const & stateMasks, vector<unsigned> & maxVariables)
  {
    // init data structures
    if (inMessages_.size() == 0)
      initMessages();
    if (maxNeighborStates_.size() == 0)
      initMaxNeighbourStates();
    
    return runMaxSum(stateMasks, maxVariables, inMessages_, outMessages_, maxNeighborStates_);
  }


  number_t DFG::runMaxSum(stateMaskVec_t const & stateMasks, vector<unsigned> & maxVariables, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const
  {
    if(roots.size() == 0)
      errorAbort("DiscreFactorGraph.cpp::runMaxSum: No root nodes");
    
    if (maxVariables.size() == 0)
      initMaxVariables(maxVariables);

    number_t maxVal = 1;

    for(int i = 0; i < roots.size(); ++i){
      unsigned root = roots.at(i);
      runMaxSumInwardsRec(root, root, stateMasks, inMessages, outMessages, maxNeighborStates);

      // find max at root node
      unsigned dim = nodes[root].getDimension();
      message_t v(vector_t(dim), 0);
      calcSumProductMessageVariable(root, root, stateMasks[root], inMessages[root], v);  // sumProduct and maxSum perform the same calculations for variables
      unsigned maxState = ublas::index_norm_inf(v.first);
      maxVal *= v.first[maxState] * std::exp(v.second);

      //  backtrack
      backtrackMaxSumOutwardsRec(maxVariables, root, root, maxState, maxNeighborStates);
    }
    return maxVal;
  }


  void DFG::runMaxSumInwardsRec(unsigned current, unsigned sender, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    // recursively call all nodes
    for (unsigned i = 0; i < nbs.size(); i++) {
      unsigned nb = nbs[i];
      if (nb != sender) 
	runMaxSumInwardsRec(nb, current, stateMasks, inMessages, outMessages, maxNeighborStates);
    }

    // calc outMessage and store in outMessages -- inMessages provides links to outMessages, which provides the message passing.
    if (current == sender) // this is root node of recursion
      return;
    calcMaxSumMessage(current, sender, stateMasks, inMessages, outMessages, maxNeighborStates); //  receiver = sender
  }


  void DFG::calcMaxSumMessage(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const
  {
    if (nodes[current].isFactor())
      calcMaxSumMessageFactor(current, receiver, inMessages, outMessages, maxNeighborStates);
    else
      calcSumProductMessageVariable(current, receiver, stateMasks, inMessages, outMessages); // calc for variables same as in sum product algorithm
  };


  void DFG::calcMaxSumMessageFactor(unsigned current, unsigned receiver, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages, vector<vector<vector<unsigned> > > & maxNeighborStates) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    message_t & outMes = outMessages[current][ getIndex( nbs, receiver) ];  // identify message
    vector< message_t const *> const & inMes( inMessages[current] );
    vector<vector<unsigned> > & maxNBStates = maxNeighborStates[ convNodeToFac(current) ];

    calcMaxSumMessageFactor(current, receiver, inMes, outMes, maxNBStates);
  }


  void DFG::calcMaxSumMessageFactor(unsigned current, unsigned receiver, vector<message_t const *> const & inMes, message_t & outMes, vector<vector<unsigned> > & maxNBStates) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    DFGNode const & nd = nodes[current];
    matrix_t const & pot = nd.getPotential();

    // one neighbor 
    if (nd.getDimension() == 1) {
      for (unsigned i = 0; i < pot.size2(); i++)
	outMes.first[i] = pot(0, i);
      return;
    }

    // two neighbors
    if (nd.getDimension() == 2) {
      if (nbs[0] == receiver) {
	unsigned nbIdx = 1;
	for (unsigned i = 0; i < pot.size1(); i++) {
	  message_t v(elemProd<vector_t>( inMes[nbIdx]->first, ublas::matrix_row<matrix_t const>(pot, i)), inMes[nbIdx]->second);
	  unsigned maxIdx = ublas::index_norm_inf(v.first);
	  maxNBStates[nbIdx][i] = maxIdx;
	  outMes.first[i] = v.first[maxIdx];
	  outMes.second = v.second;
	}
      }
      else { // nbs[1] == receiver
	unsigned nbIdx = 0;
	for (unsigned i = 0; i < pot.size2(); i++) {
	  message_t v( elemProd<vector_t>( inMes[nbIdx]->first, ublas::matrix_column<matrix_t const>(pot, i) ), inMes[nbIdx]->second);
	  unsigned maxIdx = ublas::index_norm_inf(v.first);
	  maxNBStates[nbIdx][i] = maxIdx;
	  outMes.first[i] = v.first[maxIdx];
	  outMes.second = v.second;
	}
      }
      return;
    }
    
    // more than two neighbors -- should not happen
    errorAbort("calcMaxSumMessageFactor: Error, factor with more than two neighbors");
  }


  void DFG::backtrackMaxSumOutwardsRec(vector<unsigned> & maxVariables, unsigned current, unsigned sender, unsigned maxState, vector<vector<vector<unsigned> > > & maxNeighborStates) const
  {
    // store variable max state
    if (not nodes[current].isFactor())
      maxVariables[ convNodeToVar(current) ] = maxState;

    vector<unsigned> const & nbs = neighbors[current];
    // recursively call all nodes
    for (unsigned i = 0; i < nbs.size(); i++) {
      unsigned nb = nbs[i];
      if (nb != sender) {
	unsigned nbMaxState;
	if (nodes[current].isFactor())
	  nbMaxState = maxNeighborStates[ convNodeToFac(current) ][i][maxState];
	else
	  nbMaxState = maxState; 
	backtrackMaxSumOutwardsRec(maxVariables, nb, current, nbMaxState, maxNeighborStates); //  receiver = nb
      }
    }
  }


  number_t DFG::calcNormConst(stateMaskVec_t const & stateMasks)
  {
    if (inMessages_.size() == 0)
      initMessages();

    return calcNormConst(stateMasks, inMessages_, outMessages_);
  }


  number_t DFG::calcNormConst(stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages) const
  {
    number_t res = 1;
    if(roots.size() == 0)
      errorAbort("DiscreFactorGraph.cpp::calcNormConst: No root nodes");

    for(int i = 0; i < roots.size(); ++i){
      unsigned const root = roots.at(i);
      runSumProductInwardsRec(root, root, stateMasks, inMessages, outMessages);
      res *= calcNormConstComponent(root, stateMasks[root], inMessages[root]);
    }
    return res;
  }


  number_t DFG::calcNormConstComponent(unsigned varId, stateMask_t const * stateMask, vector<message_t const *> const & inMes) const
  {
    // Calculates Normalizing Constant for the component that contains varId
    assert( varId < variables.size() );
    unsigned dim = nodes[ variables[varId] ].getDimension();
    message_t v(vector_t(dim), 0);
    calcSumProductMessageVariable(varId, varId, stateMask, inMes, v);
    return ublas::sum(v.first)*std::exp(v.second);
  }

  void DFG::calcVariableMarginals(stateMaskVec_t const & stateMasks){
    calcVariableMarginals( variableMarginals_, stateMasks);
  }

  void DFG::calcVariableMarginals(vector<vector_t> & variableMarginals, stateMaskVec_t const & stateMasks)
  {
    calcVariableMarginals(variableMarginals, stateMasks, inMessages_);
  }


  void DFG::calcVariableMarginals(vector<vector_t> & variableMarginals, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMessages) const
  {
    if (variableMarginals.size() == 0)
      initVariableMarginals(variableMarginals);

    for (unsigned i = 0; i < variables.size(); i++) {
      message_t v( vector_t(variableMarginals[i].size()), 0);
      calcSumProductMessageVariable(i, i, stateMasks[i], inMessages[i], v);
      variableMarginals[i] = v.first;
    }
  }

  void DFG::calcFactorMarginals(){
    calcFactorMarginals( factorMarginals_);
  }

  void DFG::calcFactorMarginals(vector<matrix_t> & factorMarginals)
  {
    calcFactorMarginals(factorMarginals, inMessages_);
  }


  void DFG::calcFactorMarginals(vector<matrix_t> & factorMarginals, vector<vector<message_t const *> > & inMessages) const
  {
    if (factorMarginals.size() == 0)
      initFactorMarginals(factorMarginals);

    for (unsigned facId = 0; facId < factors.size(); facId++) {
      unsigned ndId = convFacToNode(facId);
      vector<message_t const *> const & inMes = inMessages[ndId];
      matrix_t & m = factorMarginals[facId];
      m = nodes[ndId].getPotential();  // set facMar = potential
      if (nodes[ndId].getDimension() == 1)
	for (unsigned i = 0; i < m.size2(); i++)
	  m(0, i) *= inMes[0]->first[i];
      else if (nodes[ndId].getDimension() == 2)
	for (unsigned i = 0; i < m.size1(); i++)
	  for (unsigned j = 0; j < m.size2(); j++)
	    m(i, j) *= inMes[0]->first[i] * inMes[1]->first[j];
      else // should not happen
	errorAbort("Factor with more than two neighbors. Aborts.");
    }

    //Not same normalization constant for all factors when graph is disconnected
    for (unsigned facId = 0; facId < factors.size(); facId++){
      number_t const Z = sumMatrix(factorMarginals[facId]);
      factorMarginals[facId] *= 1/ Z;
    }
  }

  // init data structures given in DFG.
  void DFG::initMessages()
  {
    initMessages(inMessages_, outMessages_);
  }


  void DFG::initMaxNeighbourStates()
  {
    initMaxNeighbourStates(maxNeighborStates_);
  }


  void DFG::initMessages(vector<vector<message_t const *> > & inMessages, vector<vector<message_t> > & outMessages) const
  {
    // outMessages
    // outMessages[i][j] is a the outmessage for the j'th neighbor of node i. It's dimension equals the dimension of the linking variable node.
    outMessages.clear();
    outMessages.resize( nodes.size() );

    for (unsigned i = 0; i < nodes.size(); i++) {
      DFGNode const & nd = nodes[i];
      if (nd.isFactor()) {
	for (unsigned j = 0; j < neighbors[i].size(); j++) {
	  unsigned const dim = nodes[ neighbors[i][j] ].getDimension(); // all neighbors are variables
	  outMessages[i].push_back( message_t(vector_t(dim), 0) ); 
	}
      }
      else {  // is variable
	unsigned const dim = nodes[i].getDimension(); 
	for (unsigned j = 0; j < neighbors[i].size(); j++)
	  outMessages[i].push_back( message_t(vector_t(dim), 0) ); 
      }
    }

    // inMessages
    // inMessages[i][j] is the inMessage from the j'th neighbor of node i. It is a pointer to the corresponding outMessage. (The dimension of the vector pointed to equals the linking variable node).
    inMessages.clear();
    inMessages.resize( nodes.size() );
    for (unsigned i = 0; i < nodes.size(); i++)
      for (unsigned j = 0; j < neighbors[i].size(); j++) {
	unsigned nb = neighbors[i][j];
	unsigned k  = getIndex(neighbors[nb], i);
	inMessages[i].push_back(& outMessages[nb][k]);
      }
  }
  

  void DFG::initVariableMarginals(vector<vector_t> & variableMarginals) const
  {
    initGenericVariableMarginals(variableMarginals, *this);
  }

  void DFG::initFactorMarginals(vector<matrix_t> & factorMarginals) const
  {
    initGenericFactorMarginals(factorMarginals, *this);
  }

  void DFG::initComponents(){
    components.resize( nodes.size(), 0 );

    unsigned currentComponent = 0;
    for(int i = 0; i < components.size(); ++i){
      if(components.at(i) != 0)
	continue; //already visited

      //Do BFS starting at node i      
      currentComponent++;
      std::queue<unsigned> toVisit;
      toVisit.push(i);
      roots.push_back(i);

      while(!toVisit.empty()){
	//Visiting
	int visit = toVisit.front();
	toVisit.pop();
	components.at(visit) = currentComponent;

	//Push neighbors on queue
	for(int j = 0; j < neighbors.at(visit).size(); ++j){
	  //check if already visited
	  if( components.at( neighbors.at(visit).at(j) ) == 0)
	    toVisit.push( neighbors.at(visit).at(j) );
	}
      }
      

    }
  }

  // IS sample from factor graph
  void DFG::sampleIS(boost::mt19937 & gen, vector<vector_t> const & varMarginals, vector<matrix_t> const & facMarginals, vector<vector_t> const & ISVarMarginals, vector<matrix_t> const & ISFacMarginals, vector<unsigned> & sim, number_t & weight){
    //Make assertions
    if( varMarginals.size() != variables.size() )
      errorAbort("DFG::sampleIS: varMarginals.size() != variables.size()");
    if( facMarginals.size() != factors.size() )
      errorAbort("DFG::sampleIS: facMarginals.size() != factors.size()");
    if( ISVarMarginals.size() != variables.size() )
      errorAbort("DFG::sampleIS: ISVarMarginals.size() != variables.size()");
    if( ISFacMarginals.size() != factors.size() )
      errorAbort("DFG::sampleIS: ISFacMarginals.size() != factors.size()");
    
    //return vector
    sim.resize(variables.size());
    weight = 1;

    for(int r = 0; r < roots.size(); ++r){
      unsigned root = roots.at(r);
      
      if(nodes.at(root).isFactor())
	errorAbort("DiscreFactorGraph.cpp::DFG::sampleIS: Root nodes currently have to be variables");

      unsigned varId = convNodeToVar(root);
      boost::random::discrete_distribution<int, number_t> dist( ISVarMarginals.at(varId).begin(), ISVarMarginals.at(varId).end());
      boost::variate_generator<boost::mt19937 &, boost::random::discrete_distribution< int, number_t> > distGen( gen, dist);
      
      unsigned state = distGen();
      sim.at(varId) = state;
      weight *= varMarginals.at(varId)(state)/ISVarMarginals.at(varId)(state)/sum(varMarginals.at(varId))*sum(ISVarMarginals.at(varId));

      simulateVariableIS(gen, root, root, state, facMarginals, ISFacMarginals, sim, weight);
    }
  }

    void DFG::simulateVariableIS(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<matrix_t> const & ISFacMarginals, vector<unsigned> & sim, number_t & weight){
      //Loop over neighbors except sender ( notice if current==sender, this is the root node)
      vector<unsigned> const & nbs = neighbors.at(current);
      for(int i = 0; i < nbs.size(); ++i){
	unsigned nb = nbs.at(i);
	if(nb != sender)
	  simulateFactorIS(gen, nb, current, state, facMarginals, ISFacMarginals, sim, weight);
      }
    }

  void DFG::simulateFactorIS(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<matrix_t> const & ISFacMarginals,  vector<unsigned> & sim, number_t & weight){
    //At most a single nb except root
    vector<unsigned> const & nbs = neighbors.at(current);

    if( nbs.size() == 1)
      return; //Prior for some variable that has already been set
    else if( nbs.size() == 2){
      unsigned facId = convNodeToFac(current);

      if(nbs[0] == sender){
	const boost::numeric::ublas::matrix_row<const matrix_t> prob( facMarginals.at(facId), state);
	const boost::numeric::ublas::matrix_row<const matrix_t> ISprob( ISFacMarginals.at(facId), state);
	boost::random::discrete_distribution<int, number_t> dist( ISprob.begin(), ISprob.end() );
	boost::variate_generator<boost::mt19937 &, boost::random::discrete_distribution<int, number_t> > distGen( gen, dist);
	unsigned receiver_state = distGen();
	
	//set state
	sim.at( convNodeToVar(nbs[1])) = receiver_state;
	
	//set weight
	weight *= prob(receiver_state)/ISprob(receiver_state)/sum(prob)*sum(ISprob);

	//call neighbour
	simulateVariableIS(gen, nbs[1], current, receiver_state, facMarginals, ISFacMarginals, sim, weight);
      }
      else if( nbs[1] == sender){
	const boost::numeric::ublas::matrix_column<const matrix_t> prob( facMarginals.at(facId), state);
	const boost::numeric::ublas::matrix_column<const matrix_t> ISprob( ISFacMarginals.at(facId), state);
	boost::random::discrete_distribution<int, number_t> dist( ISprob.begin(), ISprob.end() );
	boost::variate_generator<boost::mt19937 &, boost::random::discrete_distribution<int, number_t> > distGen( gen, dist);
	unsigned receiver_state = distGen();
	
	//set state
	sim.at( convNodeToVar(nbs[0])) = receiver_state;

	//set weight
	weight *= prob(receiver_state)/ISprob(receiver_state)/sum(prob)*sum(ISprob);
	
	//call neighbour
	simulateVariableIS(gen, nbs[0], current, receiver_state, facMarginals, ISFacMarginals, sim, weight);
      }
    }
    else{
      errorAbort("DiscreFactorGraph.cpp::DFG::simulateFactor: Reached factor with " + toString(nbs.size()) + " neighbors, only 1 or 2 allowed");
    }    
  }

  vector<unsigned> DFG::sample(boost::mt19937 & gen, vector<vector_t> const & varMarginals, vector<matrix_t> const & facMarginals){
    vector<unsigned> ret;
    sample(gen, varMarginals, facMarginals, ret);
    return ret;
  }

  // Sample from factor graph
  void DFG::sample(boost::mt19937 & gen, vector<vector_t> const & varMarginals, vector<matrix_t> const & facMarginals, vector<unsigned> & sim){
    
    //TODO: Do some input checks!

    //Return vector
    sim.resize(variables.size());

    for(int r = 0; r < roots.size(); ++r){
      unsigned root = roots.at(r);
      
      //Sample from the root marginal distribution. Should be variable
      if(nodes.at(root).isFactor())
      	errorAbort("DiscreFactorGraph.cpp::DFG::sample: Root nodes currently have to be variables");

      unsigned varId = convNodeToVar(root);
      boost::random::discrete_distribution<int, number_t> dist( varMarginals.at(varId).begin(), varMarginals.at(varId).end());
      boost::variate_generator<boost::mt19937 &, boost::random::discrete_distribution< int, number_t> > distGen( gen, dist);

      unsigned state = distGen();
      sim.at(varId) = state;

      //Call simulateVariable
      simulateVariable(gen, root, root, state, facMarginals, sim);
    }

  }

  void DFG::simulateVariable(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<unsigned> & sim){
    //Loop over neighbors except sender (notice if current==sender, this is the root node)
    vector<unsigned> const & nbs = neighbors.at(current);
    for(int i = 0; i < nbs.size(); ++i){
      unsigned nb = nbs.at(i);
      if(nb != sender)
	simulateFactor(gen, nb, current, state, facMarginals, sim);
    }
  }

  void DFG::simulateFactor(boost::mt19937 & gen, unsigned current, unsigned sender, unsigned state, vector<matrix_t> const & facMarginals, vector<unsigned> & sim){
    //At most a single nb except root
    vector<unsigned> const & nbs = neighbors.at(current);

    if( nbs.size() == 1)
      return; //Prior for some variable that has already been set
    else if( nbs.size() == 2){
      //Two neighbors send to other
      unsigned facId = convNodeToFac(current);
      if(nbs[0] == sender){
        const boost::numeric::ublas::matrix_row<const matrix_t> prob( facMarginals.at(facId), state);
	boost::random::discrete_distribution<int, number_t> dist( prob.begin(), prob.end() );
	boost::variate_generator<boost::mt19937 &, boost::random::discrete_distribution<int, number_t> > distGen( gen, dist);
	unsigned receiver_state = distGen();
	
	//set state
	sim.at( convNodeToVar(nbs[1])) = receiver_state;
	
	//call neighbour
	simulateVariable(gen, nbs[1], current, receiver_state, facMarginals, sim);
      }
      else if( nbs[1] == sender){
	const boost::numeric::ublas::matrix_column<const matrix_t> prob( facMarginals.at(facId), state);
	boost::random::discrete_distribution<int, number_t> dist( prob.begin(), prob.end() );
	boost::variate_generator<boost::mt19937 &, boost::random::discrete_distribution<int, number_t> > distGen( gen, dist);
	unsigned receiver_state = distGen();
	
	//set state
	sim.at( convNodeToVar(nbs[0])) = receiver_state;
	
	//call neighbour
	simulateVariable(gen, nbs[0], current, receiver_state, facMarginals, sim);
      }
    }
    else{
      errorAbort("DiscreFactorGraph.cpp::DFG::simulateFactor: Reached factor with " + toString(nbs.size()) + " neighbors, only 1 or 2 allowed");
    }
  }


  
  // get max neighbor dimension
  unsigned DFG::maxNeighborDimension(vector<unsigned> const & nbs) const
  {
    unsigned maxDim = 0;
    BOOST_FOREACH(unsigned idx, nbs) {
      unsigned dim = nodes[idx].getDimension();
      if ( dim > maxDim) 
	maxDim = dim;
    }
    return maxDim;
  }


  // Sets up the data structure needed for maxSum. Indexing:
  // maxNeighborStates[i][j] returns a vector of unsigned (states)
  // corresponding to the j'th neigbor of the i'th factor.
  void DFG::initMaxNeighbourStates(vector<vector<vector<unsigned> > > & maxNeighborStates) const
  {
    maxNeighborStates.clear();

    maxNeighborStates.resize( factors.size() );
    for (unsigned i = 0; i < factors.size(); i++) {
      unsigned const dim = nodes[ factors[i] ].getDimension();
      maxNeighborStates[i].resize(dim);
      unsigned maxVarDim = maxNeighborDimension( getFactorNeighbors(i) );
      for (unsigned j = 0; j < dim; j++)
	maxNeighborStates[i][j].resize(maxVarDim);
    }
  }
  

  void DFG::initMaxVariables(vector<unsigned> & maxVariables) const
  {
    maxVariables.clear();
    maxVariables.resize( variables.size() );
  }

  vector<number_t> calcNormConsMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    vector<number_t> result( stateMask2DVec.size() );
    calcNormConsMultObs(result, stateMask2DVec, dfg);
    return result;
  }

  void DFG::init(vector<unsigned> const & varDimensions, vector<matrix_t> const & facPotentials, vector<vector<unsigned> > const & facNeighbors)
  {
    // reserve memory
    nodes.reserve( varDimensions.size() + facPotentials.size() ); 
    neighbors.reserve( varDimensions.size() + facPotentials.size() ); 
    variables.reserve( varDimensions.size() );
    factors.reserve( facPotentials.size() );

    // define variable nodes
    unsigned idx = 0;
    BOOST_FOREACH(unsigned dim, varDimensions) {
      nodes.push_back( DFGNode(dim) );
      variables.push_back(idx);
      idx++;
    }
      
    // define factor nodes
    BOOST_FOREACH(matrix_t const & pot, facPotentials) {
      nodes.push_back( DFGNode( pot ) );
      factors.push_back(idx);
      idx++;
    }

    // runSumProduct
    

    neighbors.resize( variables.size() ); // placeholder for var neighbors
    neighbors.insert(neighbors.end(), facNeighbors.begin(), facNeighbors.end());  // all facNeighbors refer to variables, which have not changed indices
    addTwoWayLinks(neighbors); // now add the var 

    initComponents();
    consistencyCheck();
  }

  //Calculate full likelihood
  number_t DFG::calcFullLikelihood( vector<unsigned> const & sample){
    //Make assertions
    if( sample.size() != variables.size())
      errorAbort("DFG::calcFullLikelihood: Not full sample");
    number_t res = 1; 

    //Loop over all factors
    for(int facId = 0; facId < factors.size(); ++facId){
      //Find neighbors      
      vector<unsigned> const & nbs = neighbors[ convFacToNode(facId) ];

      //convtovar read in sample 
      if( nbs.size() == 1)
	res *= nodes[ convFacToNode(facId) ].getPotential()(0, sample.at( convNodeToVar(nbs[0]) ));
      if( nbs.size() == 2)
	res *= nodes[ convFacToNode(facId) ].getPotential()( sample.at( convNodeToVar(nbs[0])), sample.at( convNodeToVar(nbs[1])) );
    }
    return res;
  }

  //Calculate expectancies
  pair<number_t,number_t> DFG::calcExpect(vector<matrix_t> const & fun_a, vector<matrix_t> const & fun_b, stateMaskVec_t const & stateMasks){
    //Initialize messages
    if( inMu_.size() == 0)
      initMessages(inMu_, outMu_);
    if( inLambda_.size() == 0)
      initMessages(inLambda_, outLambda_);

    number_t res_lik = 1; // Likelihood type
    number_t res_exp = 0; //Expectancy type

    for(int i = 0; i < roots.size(); ++i){
      //Calc inward recursion
      unsigned root = roots.at(i);
      runExpectInwardsRec(root, root, fun_a, fun_b, stateMasks, inMu_, outMu_, inLambda_, outLambda_);

      //use incoming messages to root to calculate expectancy
      stateMask_t const * stateMask = stateMasks[ convNodeToVar(root) ];
      
      // calculate the likelihood of current component
      number_t lik_com = calcNormConstComponent(root, stateMask, inMu_[root]);
      res_lik *= lik_com;

      if(!nodes[root].isFactor()){
	unsigned dim = nodes[root].getDimension();
	message_t mu(vector_t(dim), 0);
	message_t lambda(vector_t(dim), 0);

	calcExpectMessageVariable(root, root, stateMask, inMu_[root], mu, inLambda_[root], lambda);
	res_exp += sum(lambda.first)*std::exp(mu.second)/lik_com;
      }
    
    }//Ends loop over roots

    return make_pair(res_lik, res_exp*res_lik);
  }

  void DFG::runExpectInwardsRec(unsigned current, unsigned sender, vector<matrix_t> const & fun_a, vector<matrix_t> const & fun_b, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMu, vector<vector<message_t> > & outMu, vector<vector<message_t const *> > & inLambda, vector<vector<message_t> > & outLambda) const{
    // recursively call all nodes
    vector<unsigned> const & nbs = neighbors[current];
    for(unsigned i = 0; i < nbs.size(); ++i){
      unsigned nb = nbs[i];
      if(nb != sender)
	runExpectInwardsRec(nb, current, fun_a, fun_b, stateMasks, inMu, outMu, inLambda, outLambda);
    }

    if( current == sender) // this is root so does not send anything inwards
      return;
    
    calcExpectMessage(current, sender, fun_a, fun_b, stateMasks, inMu, outMu, inLambda, outLambda);
  }

  void DFG::calcExpectMessage(unsigned current, unsigned receiver, vector<matrix_t> const & fun_a, vector<matrix_t> const & fun_b, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMu, vector<vector<message_t> > & outMu, vector<vector<message_t const *> > & inLambda, vector<vector<message_t> > & outLambda) const {
    if(nodes[current].isFactor())
      calcExpectMessageFactor(current, receiver, fun_a, fun_b, inMu, outMu, inLambda, outLambda);
    else
      calcExpectMessageVariable(current, receiver, stateMasks, inMu, outMu, inLambda, outLambda);
  }

  void DFG::calcExpectMessageFactor(unsigned current, unsigned receiver, vector<matrix_t> const & fun_a, vector<matrix_t> const & fun_b, vector<vector<message_t const *> > & inMu, vector<vector<message_t> > & outMu, vector<vector<message_t const *> > & inLambda, vector<vector<message_t> > & outLambda) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    vector< message_t const *> const & inMesMu( inMu[current] );
    message_t & outMesMu = outMu[current][ getIndex( nbs, receiver) ]; //identify message
    vector< message_t const *> const & inMesLambda( inLambda[current] );
    message_t & outMesLambda = outLambda[current][ getIndex( nbs, receiver) ]; //identify message

    calcExpectMessageFactor(current, receiver, fun_a, fun_b, inMesMu, outMesMu, inMesLambda, outMesLambda);
  }

  void DFG::calcExpectMessageFactor(unsigned current, unsigned receiver, vector<matrix_t> const & fun_a, vector<matrix_t> const & fun_b, vector<message_t const *> const & inMesMu, message_t & outMesMu, vector<message_t const *> const & inMesLambda, message_t & outMesLambda) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    DFGNode const & nd = nodes[current];
    matrix_t pot_a = fun_a.at( convNodeToFac(current));
    matrix_t pot_b = fun_b.at( convNodeToFac(current));

    //One neighbor
    if(nd.getDimension() == 1){
      for( unsigned i = 0; i < nd.getPotential().size2(); ++i){
	outMesMu.first[i] = pot_a(0,i);
	outMesLambda.first[i] = pot_a(0,i)*pot_b(0,i);
	outMesMu.second = 0; // Use same normalization constant for both types of messages
      }
    }

    //Two neighbors
    if(nd.getDimension() == 2){
      if(nbs[0] == receiver){ // rows are receivers
	outMesMu.first = prod(pot_a, inMesMu[1]->first);
	outMesMu.second = inMesMu[1]->second;
	
	//TODO Vectorize the following computations
	for(unsigned i = 0; i < pot_a.size1(); ++i){
	  outMesLambda.first[i] = 0;
	  for(unsigned j = 0; j < pot_a.size2(); ++j){
	    outMesLambda.first[i] += pot_a(i,j) * pot_b(i,j) * inMesMu[1]->first[j];
	    outMesLambda.first[i] += pot_a(i,j) * inMesLambda[1]->first[j];
	  }
	}
      }
      else{ // nbs[1] == receiver, columns are receivers
	outMesMu.first = prod(inMesMu[0]->first, pot_a);
	outMesMu.second = inMesMu[0]->second;

	//TODO Vectorize the following computations
	for(unsigned j = 0; j < pot_a.size2(); ++j){
	  outMesLambda.first[j] = 0;
	  for(unsigned i = 0; i < pot_a.size1(); ++i){
	    outMesLambda.first[j] += pot_a(i,j) * pot_b(i,j) * inMesMu[0]->first[i];
	    outMesLambda.first[j] += pot_a(i,j) * inMesLambda[0]->first[i];
	  }
	}
      }
    }
  }

  void DFG::calcExpectMessageVariable(unsigned current, unsigned receiver, stateMaskVec_t const & stateMasks, vector<vector<message_t const *> > & inMu, vector<vector<message_t> > & outMu, vector<vector<message_t const *> > & inLambda, vector<vector<message_t> > & outLambda) const
  {
    vector<unsigned> const & nbs = neighbors[current];
    vector< message_t const *> const & inMesMu( inMu[current]); //Relies on precondition, possibly do assertion
    message_t & outMesMu = outMu[current][ getIndex( nbs, receiver) ];
    vector< message_t const *> const & inMesLambda( inLambda[current]);
    message_t & outMesLambda = outLambda[current][ getIndex( nbs, receiver) ];

    stateMask_t const * stateMask = stateMasks[ convNodeToVar(current) ];

    calcExpectMessageVariable(current, receiver, stateMask, inMesMu, outMesMu, inMesLambda, outMesLambda);
  }

  void DFG::calcExpectMessageVariable(unsigned current, unsigned receiver, stateMask_t const * stateMask, vector<message_t const *> const & inMesMu, message_t & outMesMu, vector<message_t const *> const & inMesLambda, message_t & outMesLambda) const
  {
    vector<unsigned> const & nbs = neighbors[current];

    //Calculate mu messages
    outMesMu.second = 0;
    if (stateMask){
      for(unsigned i = 0; i < outMesMu.first.size(); ++i)
	outMesMu.first[i] = (*stateMask)[i];
    }
    else{ //init to 1
      for(unsigned i = 0; i < outMesMu.first.size(); ++i)
	outMesMu.first[i] = 1;
    }

    for(unsigned i = 0; i < nbs.size(); ++i){
      if(nbs[i] != receiver){
	outMesMu.second += inMesMu[i]->second;
	outMesMu.first = elemProd<vector_t>( inMesMu[i]->first, outMesMu.first);
      }
    }

    // Normalize message
    double nc = sum(outMesMu.first);
    outMesMu.first /= nc;
    outMesMu.second += std::log(nc);

    // Calculate lambda messages
    outMesLambda.first *= 0;

    for(unsigned i = 0; i < nbs.size(); ++i){
      if(nbs[i] == receiver)
	continue;
      vector_t add = inMesLambda[i]->first;
      for(unsigned j = 0; j < nbs.size(); ++j){
	if(j==i or nbs[j] == receiver)
	  continue;
	add = elemProd<vector_t>(inMesMu[j]->first, add);
      }
      outMesLambda.first += add;
    }

    if(stateMask){//observed variable
      for(unsigned i = 0; i < outMesLambda.first.size(); ++i)
	outMesLambda.first[i] *= (*stateMask)[i];
    }

    // Normalize 
    outMesLambda.first /= nc;
  }

  void calcNormConsMultObs(vector<number_t> & result, stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    for (unsigned i = 0; i < stateMask2DVec.size(); i++)
      result[i] = dfg.calcNormConst(stateMask2DVec[i]);
  }


  // Calculates the accumulated variable marginals over all observation vectors. 
  vector<vector_t> calcVarAccMarMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    vector<vector_t> accVariableMarginals;
    initAccVariableMarginals(accVariableMarginals, dfg);
    calcVarAccMarMultObs(accVariableMarginals, stateMask2DVec, dfg);
    return accVariableMarginals;
  }


  void calcVarAccMarMultObs(vector<vector_t> & result, stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    unsigned varCount = dfg.variables.size();
    assert(result.size() == varCount);

    vector<vector_t> tmpVarMar;
    dfg.initVariableMarginals(tmpVarMar);
    for (unsigned i = 0; i < stateMask2DVec.size(); i++) {
      dfg.runSumProduct(stateMask2DVec[i]);
      dfg.calcVariableMarginals(tmpVarMar, stateMask2DVec[i]);
      for (unsigned j = 0; j < varCount; j++)
	result[j] += tmpVarMar[j]; 
    }
  }


  // Calculates the accumulated factor marginals over all observation vectors. 
  vector<matrix_t> calcFacAccMarMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    vector<matrix_t> accFactorMarginals;
    initAccFactorMarginals(accFactorMarginals, dfg);
    calcFacAccMarMultObs(accFactorMarginals, stateMask2DVec, dfg);
    return accFactorMarginals;
  }


  void calcFacAccMarMultObs(vector<matrix_t> & result, stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    unsigned facCount = dfg.factors.size();
    assert(result.size() == facCount);

    vector<matrix_t> tmpFacMar;
    dfg.initFactorMarginals(tmpFacMar);
    for (unsigned i = 0; i < stateMask2DVec.size(); i++) {
      dfg.runSumProduct(stateMask2DVec[i]);
      dfg.calcFactorMarginals(tmpFacMar);
      for (unsigned j = 0; j < facCount; j++)
	result[j] += tmpFacMar[j];
    }
  }


  void calcVarAndFacAccMarMultObs(vector<vector_t> & varResult, vector<matrix_t> & facResult, stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    unsigned varCount = dfg.variables.size();
    assert(varResult.size() == varCount);

    unsigned facCount = dfg.factors.size();
    assert(facResult.size() == facCount);

    vector<vector_t> tmpVarMar;
    dfg.initVariableMarginals(tmpVarMar);
    vector<matrix_t> tmpFacMar;
    dfg.initFactorMarginals(tmpFacMar);
    for (unsigned i = 0; i < stateMask2DVec.size(); i++) {
      dfg.runSumProduct(stateMask2DVec[i]);
      dfg.calcVariableMarginals(tmpVarMar, stateMask2DVec[i]);
      for (unsigned j = 0; j < varCount; j++)
	varResult[j] += tmpVarMar[j];
      dfg.calcFactorMarginals(tmpFacMar);
      for (unsigned j = 0; j < facCount; j++)
	facResult[j] += tmpFacMar[j];
    }
  }


  void calcMaxProbStatesMultObs(vector<number_t> & maxProbResult, vector<vector<unsigned> > & maxVarResult, stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    unsigned long stateMaskVecCount = stateMask2DVec.size();
    unsigned varCount = dfg.variables.size();

    // adjust result references if necessary
    if (maxProbResult.size() != stateMask2DVec.size() )
      maxProbResult.resize(stateMaskVecCount);

    if (maxVarResult.size() != stateMaskVecCount)
      maxVarResult.resize(stateMaskVecCount);
    assert(varCount > 0);
    if (maxVarResult[0].size() != varCount)
      for (unsigned i = 0; i < stateMaskVecCount; i++)
	maxVarResult[i].resize(varCount);
    
    // calculations
    for (unsigned i = 0; i < stateMaskVecCount; i++)
      maxProbResult[i] = dfg.runMaxSum(stateMask2DVec[i], maxVarResult[i]);
  }


  pair<vector<number_t>, vector<vector<unsigned> > > calcMaxProbStatesMultObs(stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    long stateMaskVecCount = stateMask2DVec.size();
    unsigned varCount = dfg.variables.size();
    vector<number_t> maxProbResult(stateMaskVecCount);
    vector<vector<unsigned> > maxVarResult(stateMaskVecCount, vector<unsigned>(varCount) );

    calcMaxProbStatesMultObs(maxProbResult, maxVarResult, stateMask2DVec, dfg);
    return make_pair(maxProbResult, maxVarResult);
  }

  // tmpFacMar is used as workspace. accFacMar must be of correct size (use dfg.initFactorMarginals(tmpFacMar) )
  void calcFacAccMarAndNormConst(vector<matrix_t> & accFacMar, vector<matrix_t> & tmpFacMar, number_t & normConst, stateMaskVec_t const & stateMaskVec, DFG & dfg)
  {
    dfg.runSumProduct(stateMaskVec);
    dfg.calcFactorMarginals(tmpFacMar);
    for (unsigned j = 0; j < dfg.factors.size(); j++)
      accFacMar[j] += tmpFacMar[j];
    normConst = dfg.calcNormConst(stateMaskVec);
  }


  void calcFacAccMarAndNormConstMultObs(vector<matrix_t> & accFacMar, vector_t & normConstVec, stateMask2DVec_t const & stateMask2DVec, DFG & dfg)
  {
    // check data structures
    long unsigned stateMaskVecCount = stateMask2DVec.size();
    assert( accFacMar.size() == dfg.factors.size() );
    if( normConstVec.size() != stateMaskVecCount )
      normConstVec.resize( stateMaskVecCount );

    // calculations
    vector<matrix_t> tmpFacMar;
    dfg.initFactorMarginals(tmpFacMar);
    for (unsigned i = 0; i < stateMaskVecCount; i++) {
      calcFacAccMarAndNormConst(accFacMar, tmpFacMar, normConstVec[i], stateMask2DVec[i], dfg);
    }
  }

  // helper functions
  void initAccVariableMarginals(vector<vector_t> & variableMarginals, DFG const & dfg)
  {
    initGenericVariableMarginals(variableMarginals, dfg);
  }
    
  void initAccFactorMarginals(vector<matrix_t> & factorMarginals, DFG const & dfg)
  {
    initGenericFactorMarginals(factorMarginals, dfg);
  }


} // namespace phy
