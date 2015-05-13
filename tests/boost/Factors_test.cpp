#define BOOST_TEST_DYN_LINK

// tested code ( the defines allow access to private data
// consider if test should only include interface not implementation
#define private public
#define protected public
#include "../../src/Factors.h"
#undef private
#undef protected

// boost test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using boost::unit_test::test_suite;
using namespace phy;

// Constants
#define EPS 1e-10

/**************************************************
 Helper Functions
 **************************************************/

matrix_t mkAndFillSquareMatrixIncr(unsigned n)
{
  matrix_t m(n, n);
  for (unsigned i = 0; i < n; i++) 
    for (unsigned j = 0; j < n; j++) 
      m(i, j) = i * n + j;
  return m;
}


matrix_t mkAndFillSquareMatrixUnif(unsigned n)
{
  matrix_t m(n, n);
  for (unsigned i = 0; i < n; i++) 
    for (unsigned j = 0; j < n; j++) 
      m(i, j) = 1.0 / (n * n);
  return m;
}

/**************************************************
 Tests
 **************************************************/

BOOST_AUTO_TEST_CASE(GlobalNormFactor_mkFactor_1) 
{
  GlobalNormFactor fac("noName", mkAndFillSquareMatrixIncr(4));
  matrix_t n = fac.mkFactor();

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK(n(i, j) == i * 4 + j);
}

BOOST_AUTO_TEST_CASE(GlobalNormFactor_optimizeParameters_1) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  matrix_t pseudo = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m , pseudo); // m will be overridden by optimization below
  fac.submitCounts(m);
  fac.optimizeParameters();
  matrix_t n = fac.mkFactor();

  BOOST_CHECK_CLOSE(sumMatrix(n), 1.0, EPS);
  BOOST_CHECK_CLOSE(n(0, 0), 0.0, EPS);
  BOOST_CHECK_CLOSE(n(3, 3), 0.125, EPS);
}

BOOST_AUTO_TEST_CASE(GlobalNormFactor_constructor_1) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  matrix_t pseudo = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m, pseudo);

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK(fac.m_(i, j) == i * 4 + j);
}


BOOST_AUTO_TEST_CASE(GlobalNormFactor_submitCounts_1) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m);

  // uniform count matrix
  matrix_t c(4, 4);
  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      c(i, j) = 1;

  fac.submitCounts(c);
  fac.optimizeParameters();
  matrix_t par = fac.m_;

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK_CLOSE(par(i, j), 1.0/16, EPS);
}  

BOOST_AUTO_TEST_CASE(GlobalNormFactor_submitCounts_2) 
{
  matrix_t m = mkAndFillSquareMatrixIncr(4);
  matrix_t pseudo = mkAndFillSquareMatrixIncr(4);
  GlobalNormFactor fac("noName", m, pseudo);

  // uniform count matrix
  matrix_t c(4, 4);
  reset(c);

  fac.submitCounts(c);
  fac.optimizeParameters();

  for (unsigned i = 0; i < 4; i++) 
    for (unsigned j = 0; j < 4; j++) 
      BOOST_CHECK_CLOSE(fac.m_(i, j), pseudo(i, j) / sumMatrix(pseudo), EPS);
}  

BOOST_AUTO_TEST_CASE(RowNormFactor_general_1) 
{
  matrix_t m(4, 4);
  reset(m);
  RowNormFactor fac("noName", m);

  matrix_t c = mkAndFillSquareMatrixIncr(4);
  fac.submitCounts(c);
  fac.optimizeParameters();

  matrix_t n = fac.mkFactor();
    BOOST_CHECK(n(0, 0) == 0.0);
  for (unsigned i = 0; i < 4; i++)
    BOOST_CHECK_CLOSE(sum( row(n,i) ), 1.0, EPS);
}

BOOST_AUTO_TEST_CASE(ColumnNormFactor_general_1) 
{
  matrix_t m(4, 4);
  reset(m);
  ColumnNormFactor fac("noName", m);

  matrix_t c = mkAndFillSquareMatrixIncr(4);
  fac.submitCounts(c);
  fac.optimizeParameters();

  matrix_t n = fac.mkFactor();
    BOOST_CHECK(n(0, 0) == 0.0);
    BOOST_CHECK_CLOSE(n(3, 0), 0.5, EPS);

  for (unsigned i = 0; i < 4; i++)
    BOOST_CHECK_CLOSE(sum( column(n,i) ), 1.0, EPS);
}
