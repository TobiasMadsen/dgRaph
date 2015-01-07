/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __utilsLinAlg_h
#define __utilsLinAlg_h

////////////////////////////////////////////////////////////////
// This header defines convenient functins for ublas vector and matrix
// manipulations.
////////////////////////////////////////////////////////////////

#include "PhyDef.h"
#include "utils.h"
#include <algorithm>

namespace phy {


  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::vector<T> & v, T x = 0);

  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::matrix<T> & m, T x = 0);

  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::symmetric_matrix<T, ublas::upper> & m, T x = 0);

  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::compressed_matrix<T> & m, T x = 0);

  /** reset all entries of vector */
  template <class T>
  inline void reset(vector<T> & v);

  /** Returns sum of all matrix entries. T is some number type.*/
  template <class T>
  inline T sumMatrix(ublas::matrix<T> const & m);

  /** element-wise vector multiplication. T is a vector containing a numeric type. */
  template<class V>
  inline V elemProd(V const & v, V const & u);

  /** convert matrix m of type M into vector v of type V in place */
  template<class M, class V>
  inline void matrixToVector(V & v, M const & m);

  /** convert matrix m into vector v */
  template<class M, class V>
  inline V matrixToVector(M const & m);

  /** convert vector v into the matrix m in place*/
  template<class M, class V>
  inline void vectorToMatrix(M & m, V const & v);

  /** convert vector v into the matrix m */
  template<class V, class M>
  inline M vectorToMatrix(V const & v, unsigned size1, unsigned size2);

  /** Return bool if the two ublas matrices of type M are identical. */
  template<class M>
  inline bool matrixEqual(M const & m1, M const & m2);

  /** Convert vector to diagonal matrix. T is a numeric type. */ 
  template<class T>
  ublas::banded_matrix<T> vectorToDiagonalMatrix(ublas::vector<T> const & v);


////////////////////////////////////////////////////////////////
// inline function definitions
////////////////////////////////////////////////////////////////

  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::vector<T> & v, T x)
  {
    for (unsigned i = 0; i < v.size(); i++)
	v[i] = x;
  }


  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::matrix<T> & m, T x)
  {
    for (unsigned i = 0; i < m.size1(); i++)
      for (unsigned j = 0; j < m.size2(); j++)
	m(i, j) = x;
  }


  /** Reset all matrix entries to x (default is 0).*/
  template <class T>
  inline void reset(ublas::symmetric_matrix<T, ublas::upper> & m, T x)
  {
    for (unsigned i = 0; i < m.size1(); i++)
      for (unsigned j = 0; j < m.size2(); j++)
	m(i, j) = x;
  }

  /** Reset all matrix entries to x (default is 0). */
  template <class T>
  inline void reset(ublas::compressed_matrix<T> & m, T x ){
    m.clear();
  }

  template <class T>
  inline void reset(vector<T> & v)
  {
    for (unsigned i = 0; i < v.size(); i++)
      reset(v[i]);
  }

  /** Returns sum of all matrix entries. T is some number type.*/
  template <class T>
  inline T sumMatrix(ublas::matrix<T> const & m)
  {
    T Z = 0;
    for (unsigned i = 0; i < m.size1(); i++)
      for (unsigned j = 0; j < m.size2(); j++)
	Z += m(i, j);
    return Z;
  }


  /** element-wise vector multiplication. T is a vector containing a numeric type. */
  template<class V>
  inline V elemProd(V const & v, V const & u)
  {
    assert( v.size() == u.size() );
    V res(v);
    for (unsigned i = 0; i < res.size(); i++)
      res[i] *= u[i];
    return res;
  }


  /** convert matrix m of type M into vector v of type V in place */
  template<class M, class V>
  inline void matrixToVector(V & v, M const & m)
  {
    unsigned size1 = m.size1();
    unsigned size2 = m.size2();
    assert(v.size() == size1 * size2);
    for (unsigned i = 0; i < size1; i++)
      for (unsigned j = 0; j < size2; j++)
	v[i * size2 + j] = m(i,j);
  }


  /** convert matrix m into vector v */
  template<class M, class V>
  inline V matrixToVector(M const & m)
  {
    vector_t v( m.size1() * m.size2() );
    matrixToVector(v, m);
    return v;
  }


  /** convert vector v into the matrix m in place*/
  template<class M, class V>
  inline void vectorToMatrix(M & m, V const & v)
  {
    unsigned size2 = m.size2();
    assert(v.size() == m.size1() * size2);
    for (unsigned k = 0; k < v.size(); k++) {
      unsigned i = k / size2;
      unsigned j = k % size2;
      m(i, j) = v[k];
    }
  }


  /** convert vector v into the matrix m */
  template<class V, class M>
  inline M vectorToMatrix(V const & v, unsigned size1, unsigned size2)
  {
    assert(v.size() == size1 * size2);
    M m(size1, size2);
    vectorToMatrix(m, v);
    return m;
  }


  /** Return bool if the two ublas matrices of type M are identical. */
  template<class M>
  inline bool matrixEqual(M const & m1, M const & m2)
  {
    if (m1.size1() != m2.size1() or m1.size2() != m2.size2() )
      return false;
    return std::equal( m1.begin1(), m1.end1(), m2.begin1() );
  }


  /** Convert vector to diagonal matrix. T is a numeric type. */ 
  template<class T>
  ublas::banded_matrix<T> vectorToDiagonalMatrix(ublas::vector<T> const & v)
  {
    unsigned dim = v.size();
    ublas::banded_matrix<T> D(dim, dim, 0, 0);
    for (unsigned i = 0; i < dim; i++)
      D(i, i) = v(i);
    return D;
  }



} // end namespace phy

#endif  //__utilsLinAlg_h
