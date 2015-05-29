#include "PhyDef.h"
#include <boost/shared_ptr.hpp>

namespace phy {
  // Forward declaration
  class Potential;

  // Typedef to pointer
  typedef boost::shared_ptr<Potential> PotentialPtr_t;

  /** Keep a potential. A potential can be shared between many factors */
  class Potential{
    // Store expectation counts
    // and keep potential
    // update function for potential
  public:
    Potential(matrix_t const & pot) : potential(pot), expCounts(pot.size1(), pot.size2(), 0) {};
    
    void clearCounts();
    void submitCounts(matrix_t const & counts);
    
    // Data members
    matrix_t potential;
    matrix_t expCounts;
  };

}
