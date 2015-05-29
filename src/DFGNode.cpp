#include "DFGNode.h"
#include <map>

namespace phy{
  DFGNode::DFGNode(unsigned dimension) : 
  isFactor_(false), dimension(dimension) {}


  DFGNode::DFGNode(matrix_t const & pot) :
    isFactor_(true), potential(pot)
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
    if(isFactor_)
      return potential;
  }

  void DFGNode::setPotential(matrix_t const & pot){
    potential = pot;
  }

  DFGNode & DFGNodeSet::operator[](std::size_t n){
    if(nodeMap[n] >= DFGNodes.size())
      std::cout << "DEBUG:: Too large node index" << std::endl;
    return DFGNodes.at( nodeMap[n] );
  }

  const DFGNode & DFGNodeSet::operator[](std::size_t n) const{
    if(nodeMap[n] >= DFGNodes.size())
      std::cout << "DEBUG:: Too large node index" << std::endl;
    return DFGNodes.at( nodeMap[n] );
  }

  void DFGNodeSet::addVariable(unsigned dim){
    DFGNodes.push_back( DFGNode(dim));
    nodeMap.push_back( DFGNodes.size() - 1);
  }

  void DFGNodeSet::addFactors( vector<matrix_t> const & pot, vector<unsigned> const & potMap){
    std::map<unsigned, unsigned> nodeToPot;
    for(int i = 0; i < potMap.size(); ++i){
      if( nodeToPot.count( potMap.at(i) ) == 0){ // New potential
	DFGNodes.push_back( DFGNode( pot[ potMap[i] ]));
	nodeMap.push_back( DFGNodes.size() - 1);
	potToNodeMap[ potMap[i] ] = DFGNodes.size() - 1;
      } 
      else{ // Existing potential
	nodeMap.push_back( potToNodeMap[ potMap[i] ]);
      }
    }
  }

}
