/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/

#include "Mixtures.h"
#include "utils.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>

namespace phy {

  NormalMixture::NormalMixture(vector_t const & means, vector_t const & vars, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv, bins, means.size()), means_(means), vars_(vars), dists_(means.size())
  {
    setDists();
  }

  NormalMixture::NormalMixture( std::istream & str, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv, bins, 0), means_(), vars_(), dists_(){
    //Read parameters
    ConvertIndexNumber is(minv_, maxv_, bins_);
    number_t lin_a = is.getToNumberAlpha();
    number_t lin_b = is.getToNumberBeta();
    
    getFeatureAndSkipLine(str, "MEANS:", means_);
    getFeatureAndSkipLine(str, "VARS:", vars_);
    if(means_.size() != vars_.size() )
      Rcpp::stop("NormalMixture: The two parameter vectors must have same length");

    states_ = means_.size();

    //Convert to index space
    means_ = element_div( means_ - boost::numeric::ublas::scalar_vector<number_t>(states_,lin_b-(maxv_-minv_)/2/bins_), boost::numeric::ublas::scalar_vector<number_t>(states_,lin_a) );
    vars_ = element_div(vars_, boost::numeric::ublas::scalar_vector<number_t>(states_, lin_a*lin_a));
    
    setDists();
  }

  NormalMixture::NormalMixture(int states, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv, bins, states), means_(states), vars_(states,bins*bins/1.96/1.96), dists_()
  {
    //Set reasonable default parameter values
    for(unsigned i = 0; i < states_; ++i)
      means_(i) = (double)bins_/states_*(0.5+i);

    setDists();
  }

  void NormalMixture::serialize(std::ostream & os) const {
    ConvertIndexNumber is(minv_, maxv_, bins_);
    number_t lin_a = is.getToNumberAlpha();
    number_t lin_b = is.getToNumberBeta();

    os << "DIST:\tNORMAL" << std::endl;
    os << "MEANS:\t" << ( boost::numeric::ublas::scalar_vector<number_t>(states_,lin_b) + element_prod(boost::numeric::ublas::scalar_vector<number_t>(states_,lin_a),means_) ) << endl;
    os << "VARS:\t" << element_prod(boost::numeric::ublas::scalar_vector<number_t>(states_, lin_a*lin_a),vars_) << endl << endl;
  }

  void NormalMixture::setDists(){
    for(int i = 0; i < states_; ++i){
      if( i >= dists_.size() )
	dists_.push_back( boost::math::normal(means_(i), std::sqrt(vars_(i))) );
      else
	dists_.at(i) = boost::math::normal(means_(i), std::sqrt(vars_(i)) );
    }
  }

  number_t NormalMixture::binProb(int i, int s) const {
    //TODO assert s is within boundaries
    return cdf( dists_.at(s), i+1) - cdf( dists_.at(s), i);
  }

  void NormalMixture::mkFactor(matrix_t &m) const {
    for(int i = 0; i < m.size1(); ++i){
      for(int j = 0; j < m.size2(); ++j ){
	//TODO Calculate only in actual observations
	m(i,j) = cdf( dists_.at(i), j+1) - cdf( dists_.at(i), j);
      }
    }
  }

  int NormalMixture::optimizeParameters(matrix_t & counts){
    //TODO Parameters should be optimized in mixture
    vector_t S(states_,0);
    vector_t USS(states_,0);
    vector_t n(states_,0);

    for ( matrix_t::iterator1 it1 = counts.begin1(); it1 != counts.end1(); ++it1){
      for( matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	n(it2.index1()) += (*it2);
	S(it2.index1()) += it2.index2() * (*it2);
	USS(it2.index1()) += it2.index2() * it2.index2() * (*it2);
      }
    }
    
    for(unsigned i = 0; i < states_; ++i){
      means_(i) = S(i)/n(i);
      vars_(i) = (USS(i)-S(i)*S(i)/n(i))/(n(i)-1);
    }
    setDists();

    return 1;
  }

  GammaMixture::GammaMixture(vector_t const & alphas, vector_t const & betas, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv,bins, alphas.size()),  alphas_(alphas), betas_(betas), dists_() 
  {
    //Possibly do some assertions?
    setDists();
  }

  GammaMixture::GammaMixture(std::istream & str, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv, bins, 0), alphas_(), betas_(), dists_()
  {
    //TODO scaling!
    getFeatureAndSkipLine(str, "ALPHAS:", alphas_);
    getFeatureAndSkipLine(str, "BETAS:", betas_);
    states_ = alphas_.size();
    setDists(); 
  }

  GammaMixture::GammaMixture(unsigned states, number_t const & minv, number_t const & maxv, unsigned const & bins):Mixture(minv, maxv, bins, states), alphas_(states,3), betas_(states), dists_()
  {
    //Set reasonable defaults
    for(int i = 0; i < states_; ++i){
      betas_(i) = 3*(i+1)/maxv_;
    }
    setDists();
  }

  void GammaMixture::serialize(std::ostream & os) const{
    os << "DIST:\tGAMMA" << std::endl;
    os << "ALPHAS:\t" << alphas_ << std::endl;
    os << "BETAS:\t" << betas_ << std::endl;
  }

  number_t GammaMixture::binProb(int i, int s) const {
    if( s >= dists_.size() ) return 0;
    return cdf( dists_.at(s), i+1) - cdf( dists_.at(s), i);
  }
  
  void GammaMixture::setDists(){
    for(int i = 0; i < states_; ++i){
      if( i >= dists_.size() )
	dists_.push_back( boost::math::gamma_distribution<>(alphas_(i), 1/betas_(i) ) );
      else
	dists_.at(i) = boost::math::gamma_distribution<>(alphas_(i), 1/betas_(i) );
    }
  }

  void GammaMixture::mkFactor(matrix_t &m) const {
    for(int i = 0; i < m.size1(); ++i){
      for(int j = 0; j < m.size2(); ++j ){
	//TODO Calculate only in actual observations
	m(i,j) = cdf( dists_.at(i), j+1) - cdf( dists_.at(i), j);
      }
    }
  }

  int GammaMixture::optimizeParameters(matrix_t & counts) {
    //TODO
    setDists();
    return 1;
  }

  BetaMixture::BetaMixture(vector_t const & alphas, vector_t const & betas, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv,bins, alphas.size()),  alphas_(alphas), betas_(betas), dists_() 
  {
    //Possibly do some assertions?
    setDists();
  }


  BetaMixture::BetaMixture(std::istream & str, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv, bins, 0), alphas_(), betas_(), dists_(){
    //TODO maybe assert minv=0 and maxv=1
    string tag1,tag2,val1,val2;
    getTagFeatureAndSkipLine(str, tag1, val1);
    getTagFeatureAndSkipLine(str, tag2, val2);

    if(tag1 == "ALPHAS:" and tag2 == "BETAS:"){
      std::stringstream val1ss(val1);
      std::stringstream val2ss(val2);
      val1ss >> alphas_;
      val2ss >> betas_;

      if(alphas_.size() != betas_.size() )
	Rcpp::stop("BetaMixture: The two parameter vectors must have same length");
    }
    else if(tag1 == "N:" and tag2 == "X:"){
      alphas_ = vector_t(1,1);
      betas_ = vector_t(1,1);

      subscriptions_.push_back(val1);
      subscriptions_.push_back(val2);
    }
    else{
      Rcpp::stop("BetaMixture: Specify either ALPHAS and BETAS or N and X");
    }

    states_ = alphas_.size();
    setDists();
  }

  BetaMixture::BetaMixture(unsigned states, number_t const & minv, number_t const & maxv, unsigned const & bins) : Mixture(minv, maxv, bins, states), alphas_(states), betas_(states, 2), dists_(){
    //Set reasonable defaults
    for(unsigned i = 0; i < states_; ++i)
      alphas_(i) = i+1;

    setDists();
  }

  void BetaMixture::serialize(std::ostream & os) const {
    os << "DIST:\tBETA" << std::endl;
    if(subscriptions_.size() > 1){
      os << "N:\t" << subscriptions_.at(0) << std::endl;
      os << "X:\t" << subscriptions_.at(1) << std::endl << std::endl;
    }
    else{
      os << "ALPHAS:\t" << alphas_ << std::endl;
      os << "BETAS:\t" << betas_ << std::endl << std::endl;
    }
  }

  void BetaMixture::setDists(){
    for(int i = 0; i < states_; ++i){
      if( i >= dists_.size() )
	dists_.push_back( boost::math::beta_distribution<>(alphas_(i), betas_(i) ));
      else
	dists_.at(i) = boost::math::beta_distribution<>(alphas_(i), betas_(i));
    }
  }

  void BetaMixture::mkFactor(matrix_t &m) const {
    for(int i = 0; i < m.size1(); ++i){
      for(int j = 0; j < m.size2(); ++j){
	//We do not assume that minv=0 and maxv=1
	m(i,j) = cdf( dists_.at(i), (maxv_-minv_)*(j+1)/bins_+minv_) - cdf( dists_.at(i), (maxv_-minv_)*j/bins_+minv_);
      }
    }
  }

  int BetaMixture::optimizeParameters(matrix_t & counts){
    //Method of moments
    //TODO default to MLE
    if(getSubscriptions().size() > 0)
      return 1;

    vector_t S(states_,0);
    vector_t USS(states_,0);
    vector_t n(states_,0);
    
    for ( matrix_t::iterator1 it1 = counts.begin1(); it1 != counts.end1(); ++it1){
      for( matrix_t::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
	n(it2.index1()) += (*it2);
	S(it2.index1()) += it2.index2() * (*it2);
	USS(it2.index1()) += it2.index2() * it2.index2() * (*it2);
      }
    }

    for(unsigned i = 0; i < states_; ++i){
      //Mean
      number_t x = S(i)/n(i)/(bins_-1); //Last factor corrects for scaling
      //Variance
      number_t v = (USS(i)-S(i)*S(i)/n(i))/n(i)/(bins_-1)/(bins_-1);
      
      alphas_(i) = x*(x*(1-x)/v-1);
      betas_(i) = alphas_(i)*(1/x-1);
    }
    return 1;
  }

  void BetaMixture::update(vector<symbol_t>& var){
    //For now only implement posterior from binomial type update
    //If p is Beta(alpha, beta) distributed
    //And we observe X ~ Binomial(N, p)
    //the post. dist. of p ~ Beta(alpha+X, beta+N-X)
    //if(alphas_.size() > 1)
    //      Rcpp::Rcout << "BetaMixture::update Warning:: only first component is updated" << std::endl;
    if(var.size() > 1){
      int N = boost::lexical_cast<int>(var.at(0));
      int X = boost::lexical_cast<int>(var.at(1));
      alphas_(0) = 1 + X;
      betas_(0) = 1 + N - X;
    }
    setDists();
  }

  //IO functions
  MixPtr_t readMixture(std::istream & str, string dist, number_t const & minv, number_t const & maxv, unsigned const & bins, unsigned const & states){
    if(dist == "NORMAL" or dist == "NORM"){
      if(moreTags(str))
	return MixPtr_t( new NormalMixture(str, minv, maxv, bins) );
      else
	//Use defaults
	return MixPtr_t( new NormalMixture(states, minv, maxv, bins) );
    }
    else if(dist == "GAMMA"){
      if(moreTags(str))
	return MixPtr_t( new GammaMixture(str, minv, maxv, bins) );
      else
	//Use defaults
	return MixPtr_t( new GammaMixture(states, minv, maxv, bins) );
    }
    else if(dist == "BETA"){
      if(moreTags(str))
	return MixPtr_t( new BetaMixture(str, minv, maxv, bins) );
      else
	//Use defaults
	return MixPtr_t( new BetaMixture(states, minv, maxv, bins));
    }
    else{
      Rcpp::stop("readMixture: Unrecognized distribution type: " + dist + "\n");
    }
  }

} // end namespace phy
