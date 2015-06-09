/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Definitions_h
#define __Definitions_h

#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <utility>

namespace dgRaph {

  namespace ublas = boost::numeric::ublas;
  using namespace std;

  /** Standard float type */
  typedef double number_t ;

  /** Vector types holding number_t and xnumber_t types */
  typedef ublas::vector<number_t> vector_t;

  /** Message type, a vector with a log transformed normalization constant */
  typedef std::pair<vector_t, double> message_t;
  
  /** Matrix types holding number_t and xnumber_t types */
  typedef ublas::matrix<number_t> matrix_t;
  
  /** Type used for observed data */
  typedef std::string symbol_t;

  /** Type used to represent states of random variables */
  typedef unsigned state_t;

  /** Convert between boost::numeric::ublas::vector and std::vector */
  template<class T>
  std::vector<T> toStdVector(ublas::vector<T> v) {std::vector<T> u(v.size()); for (unsigned i = 0; i < v.size(); i++) u[i] = v[i]; return u;}

  template<class T>
  ublas::vector<T> toNumVector(std::vector<T> v) {ublas::vector<T> u(v.size()); for (unsigned i = 0; i < v.size(); i++) u[i] = v[i]; return u;}
} // end namespace dgRaph


#endif  //__Definitions_h
