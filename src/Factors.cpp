/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "Factors.h"
#include "utils.h"
#include <Rcpp.h>
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

  /** Constructor without parameter values, set to reasonable defaults */
  DiscContFactor::DiscContFactor(string const & name, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins ) : AbstractBaseFactor("discCont", name, states, bins), minv_(minv), maxv_(maxv), bins_(bins),states_(states) { 
    mixDist_ = MixPtr_t(new NormalMixture(states, minv, maxv, bins));
  }
  
  /** Constructor taking Mixture object */
  DiscContFactor::DiscContFactor(string const & name, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins, MixPtr_t mixDist ) : AbstractBaseFactor("discCont", name, states, bins), minv_(minv), maxv_(maxv), bins_(bins), states_(states), mixDist_(mixDist) {

  }

  void DiscContFactor::serialize(ostream & os) const
  {
    os << "BINS:\t" << bins_ << endl;
    os << "MIN:\t" << minv_ << endl;
    os << "MAX:\t" << maxv_ << endl;
    os << "STATES:\t" << states_ << endl;
    mixDist_->serialize(os);
  }

  int DiscContFactor::optimizeParametersImpl()
  {
    mixDist_->optimizeParameters(counts_);
  }

  void DiscContFactor::mkFactor(matrix_t &m) const
  {
    mixDist_->mkFactor(m);
  }

  vector<string> DiscContFactor::getSubscriptions() const
  {
    return mixDist_->getSubscriptions();
  }

  void ContContFactor::serialize(ostream & os) const {
    ConvertIndexNumber is1(minv1_, maxv1_, bins1_);
    ConvertIndexNumber is2(minv2_, maxv2_, bins2_);

    os << "BINS1:\t" << bins1_ << endl;
    os << "BINS2:\t" << bins2_ << endl;
    os << "MIN1:\t" << minv1_ << endl;
    os << "MAX1:\t" << maxv1_ << endl;
    os << "MIN2:\t" << minv2_ << endl;
    os << "MAX2:\t" << maxv2_ << endl;
    os << "ALPHA:\t" << alpha_* is2.getToNumberAlpha()/is1.getToNumberAlpha() << endl;
    os << "BETA:\t" << beta_*is2.getToNumberAlpha()+is2.getToNumberBeta()-alpha_*is2.getToNumberAlpha()*is1.getToNumberBeta()/is1.getToNumberAlpha() << endl;
    os << "VAR:\t" << var_*is2.getToNumberAlpha()*is2.getToNumberAlpha() << endl << endl;
  }

  void ContContFactor::mkFactor(matrix_t &m) const{
    for(matrix_t::iterator1 it1 = m.begin1(); it1 != m.end1(); ++it1){
      boost::math::normal norm( it1.index1()*alpha_+beta_, sqrt(var_) );   
      for(matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	//TODO Lazy evaluation?
	*it2 = cdf(norm, it2.index2()) - cdf(norm, it2.index2()-1);
      }
    }
  }

  int ContContFactor::optimizeParametersImpl() {
    number_t SP_xy = 0, S_x = 0, S_y = 0, USS_x = 0, USS_y = 0, N = 0;
    for(matrix_t::iterator1 it1 = counts_.begin1(); it1 != counts_.end1(); ++it1){
      for(matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	SP_xy += it2.index1()*it2.index2()  * (*it2);
	S_x += it2.index1() * (*it2);
	S_y += it2.index2() * (*it2);
	USS_x += it2.index1()*it2.index1() * (*it2);
	USS_y += it2.index2()*it2.index2() * (*it2);
	N += *it2;
      }
    }
    
    number_t SSD_x = USS_x-S_x*S_x/N;
    number_t SSD_y = USS_y-S_y*S_y/N;
    number_t SPD_xy = SP_xy - S_x*S_y/N;
    alpha_ = SPD_xy/SSD_x;
    beta_ = (S_y-S_x*alpha_)/N;
    var_ = (SSD_y - SPD_xy*SPD_xy/SSD_x)/N;
    if(var_ <= 0)
      errorAbort("ContContFactor::optimizeParametersImpl: variance 0 or less");
    
    return 1;
  }

  void BinomialFactor::serialize(ostream & os) const {
    os << "MIN:\t" << minv_ << std::endl;
    os << "MAX:\t" << maxv_ << std::endl;
    os << "PROB:\t" << prob_ << std::endl;
    os << "N:\t" << N_ << std::endl;
  }

  void BinomialFactor::mkFactor(matrix_t &m) const{
    boost::math::binomial binom(N_, prob_);
    //Return vector with probabilities over x
    for(int j = 0; j < m.size2(); ++j){
      if( j+minv_ > N_)
	m(0,j) = 0;
      else
	m(0,j) = pdf( binom, j+minv_);
    }
  }

  int BinomialFactor::optimizeParametersImpl(){
    double n = 0;
    double sum = 0;
    for(matrix_t::iterator1 it1 = counts_.begin1(); it1 != counts_.end1(); ++it1){
      for(matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	n += (*it2);
	sum += it2.index2() * (*it2);
      }
    }    
    sum += n*minv_;
    prob_ = sum/n/N_;
    return 1;
  }

  void BinomialFactor::update(vector<symbol_t>& var){
    if(updateType_ == 1 && var.size() > 0){
      //Update N_
      try{
	N_ = boost::lexical_cast<int>(var.at(0));
      }
      catch(boost::bad_lexical_cast &){
	Rcpp::Rcout << "BinomialFactor::update Warning: Variable '" << var.at(0) << "' could not be converted to int" << std::endl;
	Rcpp::stop("Failed conversion");
      }
      return;
    }
    if(updateType_ == 2 && var.size() > 0){
      //Update prob_
      try{
	prob_ = boost::lexical_cast<double>(var.at(0));
      }
      catch(boost::bad_lexical_cast &){
	Rcpp::Rcout << "BinomialFactor::update Warning: Variable '" << var.at(0) << "'could not be converted to double" << std::endl;
	Rcpp::stop("Failed conversion");
      }
      return;
    }
    if(updateType_ == 12 && var.size() > 1){
      //Update both N_ and prob_ in that order
      try{
	N_ = boost::lexical_cast<int>(var.at(0));
	prob_ = boost::lexical_cast<double>(var.at(1));
      }
      catch(boost::bad_lexical_cast &){
	Rcpp::Rcout << "BinomialFactor::update Warning: Variables '" << var.at(0) << "' and '" << var.at(1) << "'could not be converted to int and double respectively" << std::endl;
      }
      return;
    }
    Rcpp::Rcout << "BinomialFactor::update Warning: no update performed" << std::endl;
  }

  NormalMeanPostFactor::NormalMeanPostFactor(string const & name, number_t const & var, number_t const & minv, number_t const & maxv, unsigned bins, vector<string> subs) : AbstractBaseFactor("normalMeanPost", name, 1, bins), var_(var), minv_(minv), maxv_(maxv), bins_(bins)
  { 
    subscriptions = subs;
    postvar_ = 1;
    postmean_ = 0;
  }

  void NormalMeanPostFactor::serialize(ostream & os) const{
    //write out subscription factors to reflect input
    os << "BINS:\t" << bins_ << std::endl;
    os << "MIN:\t" << minv_ << std::endl;
    os << "MAX:\t" << maxv_ << std::endl;
    os << "VAR:\t" << var_ << std::endl;
    os << "X:\t";
    copy(subscriptions.begin(), subscriptions.end(), ostream_iterator<symbol_t>(os, "\t") );
    os << std::endl << std::endl;
  }

  void NormalMeanPostFactor::mkFactor(matrix_t &m) const{
    //Remember to convert to index space
    boost::math::normal dist( (postmean_-minv_)*bins_/(maxv_-minv_), 
			      std::sqrt(postvar_)*bins_/(maxv_-minv_) );
    for(int i = 0; i < m.size2(); ++i){
      m(0,i) = cdf(dist, i+1)-cdf(dist,i);
    }
  }

  void NormalMeanPostFactor::update(vector<symbol_t> & variables){
    //convert to double
    unsigned k = variables.size();
    number_t sum = 0;
    for(vector<symbol_t>::iterator it = variables.begin(); it != variables.end(); ++it){
      sum += boost::lexical_cast<double>( *it );
    }

    //find posterior mean and variance
    postmean_ = sum / k;
    postvar_ = var_ / k;
  }

  int NormalMeanPostFactor::optimizeParametersImpl(){
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
