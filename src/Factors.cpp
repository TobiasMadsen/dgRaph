/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "Factors.h"
#include "utils.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>

namespace phy {

  AbstractFullyParameterizedFactor::AbstractFullyParameterizedFactor(string const & type, string const & id, matrix_t const & m, matrix_t const & pseudoCounts) 
    : AbstractBaseFactor(type, id, m.size1(), m.size2()), m_(m), pseudoCounts_(pseudoCounts)
  {
    if (pseudoCounts.size1() != 0) {
      assert(pseudoCounts.size1() == size1_); 
      assert(pseudoCounts.size2() == size2_); 
    }
  }

  AbstractFullyParameterizedFactor::AbstractFullyParameterizedFactor(string const & type, string const & id, matrix_t const & m, matrix_t const & pseudoCounts, matrix_t const & funBMat) 
    : AbstractBaseFactor(type, id, m.size1(), m.size2()), m_(m), pseudoCounts_(pseudoCounts), funBMat_(funBMat)
  {
    if (pseudoCounts.size1() != 0) {
      assert(pseudoCounts.size1() == size1_); 
      assert(pseudoCounts.size2() == size2_); 
    }
    if (funBMat.size1() != 0) {
      assert(funBMat.size1() == size1_); 
      assert(funBMat.size2() == size2_); 
    }
  }

  void AbstractFullyParameterizedFactor::serialize(ostream& os) const{
    os << "POT_MAT:\t" << m_ << endl;
    if ( pseudoCounts_.size1() != 0 )
      os << "PC_MAT:\t"  << pseudoCounts_ << endl;
  }


  int GlobalNormFactor::optimizeParametersImpl()
  {
    m_ = counts_;
    if (pseudoCounts_.size1() != 0)
      m_ += pseudoCounts_;
    number_t n = sumMatrix(m_);
    if (n == 0.0)
      errorAbort("GlobalNormFactor:optimizeParameters: no counts: optimization not possible.");
    m_ *= 1.0 / n;
    return 1;
  }


  int ColumnNormFactor::optimizeParametersImpl()
  {
    m_ = counts_;
    if (pseudoCounts_.size1() != 0)
      m_ += pseudoCounts_;
    for (unsigned i = 0; i < size2_; i++) {
      // normalize each column
      number_t n = sum( ublas::column(m_, i) );
      if (n == 0.0)
	errorAbort("GlobalNormFactor:optimizeParameters: no counts: optimization not possible for column i=" + toString(i) + ".");
      column(m_, i) *= 1.0 / n;
    }      
    return 1;
  }
	

  int RowNormFactor::optimizeParametersImpl()
  {
    m_ = counts_;
    if (pseudoCounts_.size1() != 0)
      m_ += pseudoCounts_;
    for (unsigned i = 0; i < size1_; i++) {
      // normalize each row
      number_t n = sum( row(m_, i) );
      if (n == 0.0)
	errorAbort("GlobalNormFactor:optimizeParameters: no counts: optimization not possible for row i=" + toString(i) + ".");
      row(m_, i) *= 1.0 / n;
    }      
    return 1;
  }

  // return matrix defining factor idx 
  vector<matrix_t> AbstractBaseFactorSet::mkFactorVec() const
  {
    vector<matrix_t> v; 
    for (unsigned i = 0; i < facCount_; i++) 
      v.push_back(mkFactor(i) );
    return v;
  } 


  // set matrix defining factor idx 
  void AbstractBaseFactorSet::mkFactorVec(vector<matrix_t> & v) const 
  {
    for (unsigned i = 0; i < facCount_; i++) 
      mkFactor(v[i], i);
  }

  vector<matrix_t> AbstractBaseFactorSet::mkFunAVec() const
  {
    vector<matrix_t> v;
    for(unsigned i = 0; i < facCount_; ++i){
      v.push_back(mkFunA(i));
    }
    return v;
  }

  vector<matrix_t> AbstractBaseFactorSet::mkFunBVec() const
  {
    vector<matrix_t> v;
    for(unsigned i = 0; i < facCount_; ++i){
      v.push_back(mkFunB(i));
    }
    return v;
  }

} // end namespace phy
