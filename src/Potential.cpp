#include "Potential.h"

namespace phy {
  void Potential::clearCounts(){
    for(int i = 0; i < expCounts.size1(); ++i)
      for(int j = 0; j < expCounts.size2(); ++j)
	expCounts(i,j) = 0;
  }

  void Potential::submitCounts(matrix_t const & counts){
    expCounts += counts;
  }

}
