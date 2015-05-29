#include "PhyDef.h"
#include <vector>
#include <map>

namespace phy {

  class DFGNode
  {
  public:
    /**Variable node constructor*/
    DFGNode(unsigned dimension);
    
    /**Factor node constructor*/
    DFGNode(matrix_t const & pot);

    bool isFactor() const;

    unsigned getDimension() const;

    matrix_t getPotential() const;

    void setPotential(matrix_t const & pot);

  private:
    bool isFactor_;      // true if factor node, false if variable node
    unsigned dimension; // dimension of variable or dimension of potential
    matrix_t potential;
  };

  class DFGNodeSet{
  public:
    void addVariable(unsigned dim);
    void addFactors( vector<matrix_t> const & pot, vector<unsigned> const & potMap);

    DFGNode & operator[](std::size_t n);
    DFGNode & at(std::size_t n){ return operator [](n);}
    const DFGNode & operator[](std::size_t n) const;
    const DFGNode & at(std::size_t n) const { return operator[](n); }
    unsigned size() const { return nodeMap.size(); }

  protected:
    std::vector<DFGNode> DFGNodes;
    std::vector<unsigned> nodeMap;
    std::map<unsigned, unsigned> potToNodeMap;
  };
}
