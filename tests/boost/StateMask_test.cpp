#define BOOST_TEST_DYN_LINK

#include "../../src/StateMask.h"

// bost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using boost::unit_test::test_suite;
using namespace dgRaph;

/**************************************************
 Tests
 **************************************************/

BOOST_AUTO_TEST_CASE(StateMaskUnobserved_1){
  StateMask sm = StateMask();
  
  BOOST_CHECK_EQUAL( sm[2], 1);
}

BOOST_AUTO_TEST_CASE(StateMaskObserved_1){
  StateMaskObserved sm = StateMaskObserved(2);
  
  BOOST_CHECK_EQUAL( sm[2], 1);
  BOOST_CHECK_EQUAL( sm[1], 0);
}

BOOST_AUTO_TEST_CASE(StateMaskPosterior_1){
  vector_t post(2);
  post(0) = 0.8;
  post(1) = 0.2;

  StateMaskPosterior sm(post);
  BOOST_CHECK_EQUAL( sm[0], 0.8);
  BOOST_CHECK_EQUAL( sm[1], 0.2);
}

BOOST_AUTO_TEST_CASE(StateMaskVec_1){
  stateMaskVec_t smv(2, stateMaskPtr_t( new StateMask()) );
  smv.push_back( stateMaskPtr_t(new StateMaskObserved(2) ));
  
  BOOST_CHECK_EQUAL( (*smv[2])[1], 0);
  BOOST_CHECK_EQUAL( (*smv[2])[2], 1);
  BOOST_CHECK_EQUAL( (*smv[1])[1], 1);
  BOOST_CHECK_EQUAL( (*smv[1])[2], 1);

}
