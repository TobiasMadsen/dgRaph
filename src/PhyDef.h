/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __phyDef_h
#define __phyDef_h

#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace phy {

  namespace ublas = boost::numeric::ublas;
  using namespace std;

  /** Standard float type */
  typedef double number_t ;

  /** Float type with (potentially) extended exponent range */
  typedef long double xnumber_t ;


  /** Vector types holding number_t and xnumber_t types */
  typedef ublas::vector<number_t> vector_t;
  typedef ublas::vector<xnumber_t> xvector_t;
  
  /** Matrix types holding number_t and xnumber_t types */
  typedef ublas::matrix<number_t> matrix_t;
  typedef ublas::matrix<xnumber_t> xmatrix_t;
  
  /** Type used for observed data */
  typedef std::string symbol_t;

  /** Type used to represent states of random variables */
  typedef unsigned state_t;

  /** Type used for observed random variables */
  typedef ublas::vector<bool> stateMask_t;

  /** Type used for an enumerated set of stateMasks (e.g., input to factor graph) */
  typedef std::vector<stateMask_t const *> stateMaskVec_t;

  /** Type used for sequential data sets translated to stateMasks */
  typedef std::vector<stateMaskVec_t> stateMask2DVec_t;

  // conversion functions involving above defined types
  /** Convert extended-range-type (xnumber_t) to normal float type (number_t) */
  inline void toNumber(number_t & r, xnumber_t const & a) {r = a;}
  inline number_t toNumber(xnumber_t const & a) {return a;}

  inline void toNumber(vector_t & r, xvector_t const & v) {for (unsigned i = 0; i < v.size(); i++) r[i] = toNumber(v[i]);}
  inline vector_t toNumber(xvector_t const & v) {vector_t r(v.size()); toNumber(r, v); return r;}

  inline void toNumber(matrix_t & r, xmatrix_t const & m) {for (unsigned i = 0; i < m.size1(); i++) for (unsigned j = 0; j < m.size2(); j++) r(i, j) = toNumber( m(i, j) );}
  inline matrix_t toNumber(xmatrix_t const & m) {matrix_t r(m.size1(), m.size2()); toNumber(r, m); return r;}

  inline void toXNumber(xvector_t & r, vector_t const & v) {for (unsigned i = 0; i < v.size(); i++) r[i] = v[i];}
  inline xvector_t toXNumber(vector_t const & v) {xvector_t r(v.size()); toXNumber(r, v); return r;}

  inline void toXNumber(xmatrix_t & r, matrix_t const & m) {for (unsigned i = 0; i < m.size1(); i++) for (unsigned j = 0; j < m.size2(); j++) r(i, j) = m(i, j);}
  inline xmatrix_t toXNumber(matrix_t const & m) {xmatrix_t r(m.size1(), m.size2()); toXNumber(r, m); return r;}

  /** Convert between boost::numeric::ublas::vector and std::vector */
  template<class T>
  std::vector<T> toStdVector(ublas::vector<T> v) {std::vector<T> u(v.size()); for (unsigned i = 0; i < v.size(); i++) u[i] = v[i]; return u;}

  template<class T>
  ublas::vector<T> toNumVector(std::vector<T> v) {ublas::vector<T> u(v.size()); for (unsigned i = 0; i < v.size(); i++) u[i] = v[i]; return u;}
} // end namespace phy


#endif  //__phyDef_h
