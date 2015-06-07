#ifndef __StateMask_h
#define __StateMask_h

#include "PhyDef.h"
#include "utils.h"
#include <vector>
#include <boost/shared_ptr.hpp>

namespace phy {
  class StateMask {
  public:
    virtual const double operator[](unsigned n) const{
      return 1;
    }
  };

  class StateMaskObserved : public StateMask {
  public:
    StateMaskObserved(unsigned obs) : obs_(obs) {}

    virtual const double operator[](unsigned n) const{
      if(n == obs_)
	return 1;
      return 0;
    }

  private:
    unsigned obs_;
  };

  class StateMaskPosterior : public StateMask {
  public:
    StateMaskPosterior(vector_t const & posterior) : posterior_(posterior) {}

    virtual const double operator[](unsigned n) const{
      if(n >= posterior_.size())
	errorAbort("Wrong variable sizes. Probably specified in dataList");
      return posterior_(n);
    }

  private:
    vector_t posterior_;
  };

  /** Type used for an enumerated set of stateMasks (e.g., input to factor graph) */
  typedef boost::shared_ptr<const StateMask> stateMaskPtr_t;
  typedef std::vector<stateMaskPtr_t> stateMaskVec_t;

}

#endif // __StateMask_h
