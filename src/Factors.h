/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Factors_h
#define __Factors_h

#include "PhyDef.h"
#include "utils.h"
#include "utilsLinAlg.h"
#include "Mixtures.h"
#include <Rcpp.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

/** Classes for specifying and optimizing discrete factor potential
    (matrices). Potential matrices are specified in terms of
    parameters. Complete implementation for fully parameterized,
    normalized potentials are given. All classes allow sets of
    both variable and fixed parameters to be specified. */

/** To do: In some cases it would be convenient to be able to fix a
    subset of the otherwise free parameters in the optimization. This
    could be specified by a vector of bools, which the
    optimizationImpl could query when performing the
    optimization. */

namespace phy {

  using namespace std;

  class AbstractBaseFactor; // forward declaration

  /** Pointer to base factors */
  typedef boost::shared_ptr<AbstractBaseFactor> AbsBasFacPtr_t;

  class AbstractBaseFactor {

  public:

    virtual ~AbstractBaseFactor() {};

    /** Constructor for factor matrices. Parameters are defined in derived classes.*/
    AbstractBaseFactor(string const & type, string const & name, unsigned factorMatrixSize1, unsigned factorMatrixSize2) 
      : size1_(factorMatrixSize1), 
	size2_(factorMatrixSize2), 
	hasChanged_(false),
	type_(type),
	name_(name)
    {initCounts();}
    
    /** submit expectations for matrix entries of each factor. */
    void submitCounts(matrix_t const & counts) {counts_ += counts; hasChanged_ = true;}
										
    /** return zero 0 on failure and >0 on success. Success values may
	optionally report on other aspects of optimization. NOTE:
	derived classes should implement optimizeParametersImpl and
	NOT this function. */
    int optimizeParameters() {if (hasChanged_) {hasChanged_ = false; return optimizeParametersImpl();} else return 2;}

    /** return or set matrix defining factor */
    matrix_t mkFactor() const {matrix_t m(size1_, size2_); mkFactor(m); return m;}
    virtual void mkFactor(matrix_t & m) const = 0;

    virtual matrix_t mkFunA() const {matrix_t m(0,0); return m;}
    virtual matrix_t mkFunB() const {matrix_t m(0,0); return m;}

    /** clear all submitted expectation counts */
    void clearCounts() {reset(counts_);}

    /** returns factor type */
    string const & type() {return type_;}

    /** returns factor name */
    string const & name() const {return name_;}

    /** serialization method */
    virtual void serialize(ostream& os) const = 0;

    /** generic update method to change parameters at each observation. The factor maintains it self an update type that specifies how to parse the input and update parameters */
    virtual void update(vector<symbol_t>& var){
      Rcpp::Rcout << "AbstractBaseFactor::update: Warning: update has not been implemented for your factor" << std::endl;
    }
    
    /** get vector of variables names that the factor subscribes to */
    virtual vector<string> getSubscriptions() const{
      return subscriptions;
    }

    /** add a variable name that the factor subscribes to */
    void addSubscription(string & s){
      subscriptions.push_back(s);
    }

    void setUpdateType(int t){
      updateType_ = t;
    }

  protected:

    /** Optimize parameters implementation. Called by optimizeParameters and must be implemented by derived classes. */
    virtual int optimizeParametersImpl() = 0;  
   
    /** intialize the count matrix to the right size and all zero entries. */
    void initCounts() {counts_.resize(size1_, size2_); clearCounts();}

    /** Defines sizes of factor matrix */
    unsigned size1_;  // rows
    unsigned size2_;  // columns 

    /** stores expectation counts */
    matrix_t counts_;

    /** counts has changed -- optimization needed. */
    bool hasChanged_;

    /** Type string */
    string const type_;

    /** Name string */
    string const name_;

    /** Variables that the factor use to update parameters at each observation */
    vector<string> subscriptions;

    /** Identify which update should be performed when subscribing to variables */
    int updateType_;
  };


  class AbstractFullyParameterizedFactor : public AbstractBaseFactor {

  public:

    virtual ~AbstractFullyParameterizedFactor() {};

    /** Constructor for fully parametrized factor matrices. One
	parameter value is given for each matrix entry in matrix
	m. Optionally, one pseudoCount value can be given for each
	parameter in the matrix pseudoCount. The pseudoCount are added
	to the submitted counts before optimization. */
    AbstractFullyParameterizedFactor(string const & type, string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t() );

    AbstractFullyParameterizedFactor(string const & type, string const & name, matrix_t const & m, matrix_t const & pseudoCounts, matrix_t const & funBMat );

    virtual void mkFactor(matrix_t & m) const {m = m_;}
    using AbstractBaseFactor::mkFactor; // bringing other mkFactor definition into this name space (hidden otherwise)

    virtual matrix_t mkFunA() const { return funAMat_;}
    virtual matrix_t mkFunB() const { return funBMat_;}

    virtual void serialize(ostream& os) const;

  protected:
    matrix_t m_;
    matrix_t const pseudoCounts_;
    matrix_t funAMat_;
    matrix_t funBMat_;
    friend   void writeAbstractFullyParameterizedFactor(ostream & str, AbsBasFacPtr_t const & factorPtr);
  };


  /** matrix normalized to sum to one */
  class GlobalNormFactor : public AbstractFullyParameterizedFactor {
  public:

    virtual ~GlobalNormFactor() {};
    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    GlobalNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractFullyParameterizedFactor("globNorm", name, m, pseudoCounts) {};

    GlobalNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts, matrix_t const & funBMat) : AbstractFullyParameterizedFactor("globNorm", name, m, pseudoCounts, funBMat) {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  /** matrix columns normalized to sum to one */
  class ColumnNormFactor : public AbstractFullyParameterizedFactor  {
  public: 

    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    ColumnNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractFullyParameterizedFactor("colNorm", name, m, pseudoCounts) {};

  ColumnNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts, matrix_t const & funBMat) : AbstractFullyParameterizedFactor("colNorm", name, m, pseudoCounts, funBMat) {};
    
    virtual ~ColumnNormFactor() {};

  protected:
    virtual int optimizeParametersImpl();  
  };


  /** matrix rows normalized to sum to one */
  class RowNormFactor : public AbstractFullyParameterizedFactor  {
  public: 

    /** Constructor for fully parametrized factor matrices. The matrix values will be stored internally as floatParameters. */
    RowNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts = matrix_t()) : AbstractFullyParameterizedFactor("rowNorm", name, m, pseudoCounts) {};

  RowNormFactor(string const & name, matrix_t const & m, matrix_t const & pseudoCounts, matrix_t const & funBMat) : AbstractFullyParameterizedFactor("rowNorm", name, m, pseudoCounts, funBMat) {};

    virtual ~RowNormFactor() {};

  protected:
    virtual int optimizeParametersImpl();  
  };

  /** @brief Discrete and Continuous variable factor. Mixture of normals */
  class DiscContFactor : public AbstractBaseFactor {
  public:
    /** Constructor with parameters*/
  DiscContFactor(string const & name, vector_t const & means, vector_t const & vars, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins ) : AbstractBaseFactor("discCont", name, states, bins), minv_(minv), maxv_(maxv), bins_(bins),states_(states), mixDist_(new NormalMixture(means, vars, minv, maxv, bins)) { }

    /** Constructor using "default" parameters */
    DiscContFactor(string const & name, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins );

    /** Constructor taking Mixture object */
    DiscContFactor(string const & name, number_t const & minv, number_t const & maxv, unsigned states, unsigned bins, MixPtr_t mixDist );

    /** Destructor */
    virtual ~DiscContFactor() {} ;

    /** Serialize object(write means and variance to stream) */
    virtual void serialize(ostream & os) const;

    /** AbstractFullyParameterizedFactor just overwrite m_. Sets matrix ideally only where there are observations*/
    virtual void mkFactor(matrix_t &m) const;

    /** Pass update to Mixture */
    virtual void update( vector<symbol_t>& var){
      mixDist_->update( var);
    }

    /** Get subscriptions from Mixture */
    virtual vector<string> getSubscriptions() const;

  protected:
    virtual int optimizeParametersImpl();

  private:
    number_t minv_; ///< Endpoint of range of observations
    number_t maxv_; ///< Endpoint of range of observations
    unsigned bins_; ///< Number of bins(binning of continuous variable)
    unsigned states_; ///< Number of states
    MixPtr_t mixDist_; 
  };

  class NormalMeanPostFactor : public AbstractBaseFactor {
  public:
    NormalMeanPostFactor(string const & name, number_t const & var, number_t const & minv, number_t const & maxv_, unsigned bins, vector<string> subscriptions);

    virtual ~NormalMeanPostFactor() {} ;

    virtual void serialize(ostream & os) const;

    virtual void mkFactor(matrix_t &m) const;

    virtual void update( vector<symbol_t> & var);

  protected:
    virtual int optimizeParametersImpl();

  private:
    number_t var_;
    number_t minv_, maxv_;
    unsigned bins_;

    number_t postvar_;//Posterior variance
    number_t postmean_;//Posterior mean
  };

  /** @brief Continuous-Continuous factor linear regression of second variable in first variable
      Notice that if first variable has a normal prior, then this is effectively a 2D normal distribution */
  class ContContFactor : public AbstractBaseFactor {
  public:
    /** Constructor with parameters */
  ContContFactor(string const & name, number_t const & alpha, number_t const & beta, number_t const & var, unsigned const & bins1, unsigned const & bins2, number_t const & minv1, number_t const & maxv1, number_t const & minv2, number_t const & maxv2 ) : AbstractBaseFactor("contCont", name, bins1, bins2), alpha_(alpha), beta_(beta), var_(var), bins1_(bins1), bins2_(bins2), minv1_(minv1), maxv1_(maxv1), minv2_(minv2), maxv2_(maxv2) {}

    /** Constructor without parameters(reasonable defaults are used */
  ContContFactor(string const & name, unsigned const & bins1, unsigned const & bins2, number_t const & minv1, number_t const & maxv1, number_t const & minv2, number_t const & maxv2 ) : AbstractBaseFactor("contCont", name, bins1, bins2), alpha_(0), beta_(bins2/2), var_(bins2*bins2/1.96/1.96), bins1_(bins1), bins2_(bins2), minv1_(minv1), maxv1_(maxv1), minv2_(minv2), maxv2_(maxv2) {}

    /** Destructor */
    virtual ~ContContFactor() {} ;

    /** Serialize object(write coefficients and variance to stream) */
    virtual void serialize(ostream & os) const;

    /** Overwrite m_. Ideally set matrix only where there are observations */
    virtual void mkFactor(matrix_t &m) const;

  protected:
    virtual int optimizeParametersImpl();

  private:
    number_t alpha_, beta_, var_; ///< Regression coeffecients and variance
    unsigned bins1_, bins2_;  ///< Number of bins in each variable
    number_t minv1_,maxv1_, minv2_,maxv2_; ///< Range of observations
  };

  class BinomialFactor : public AbstractBaseFactor {
  public:
    /** Constructor minv and maxv are the ranges of the corresponding countbased statemap */
    BinomialFactor(string const & name, number_t const & prob, unsigned const & N, unsigned const & minv, unsigned const & maxv) : AbstractBaseFactor("binomial",name, 1, maxv-minv+1), N_(N), prob_(prob), minv_(minv), maxv_(maxv) { }
    
    /** Destructor */
    virtual ~BinomialFactor() {} ;
    
    /** Serialize object(write coefficients and variance to stream) */
    virtual void serialize(ostream & os) const;

    /** Overwrite m_ */
    virtual void mkFactor(matrix_t &m) const;

    /** Update parameters with new values. Useful if parameter is directly observed but is sample specific*/
    virtual void update( vector<symbol_t>& var);

  protected:
    virtual int optimizeParametersImpl();
    
  private:
    unsigned N_;
    number_t prob_; //To used for single p mode
    unsigned minv_; //min for count statemap on x
    unsigned maxv_; //max for count statemap on x
  };

  class AbstractBaseFactorSet {

  public:

    AbstractBaseFactorSet(unsigned facCount) : facCount_(facCount) {};
    virtual ~AbstractBaseFactorSet() {};

    /** submit expectations for matrix entries of each factor. */
    void submitCounts(vector<matrix_t> const & countsVec) {for (unsigned i = 0; i < facCount_; i++) submitCounts(countsVec[i], i);}
    virtual void submitCounts(matrix_t const & counts, unsigned idx) = 0;  // idx is the factor index within the class

    /** return value may optionally report on success or other aspects of optimization. */
    virtual int optimizeParameters() = 0;

    /** return or set matrix defining factor idx  */
    virtual matrix_t mkFactor(unsigned idx) const = 0; 
    virtual void mkFactor(matrix_t & m, unsigned idx) const = 0;
    virtual matrix_t mkFunA(unsigned idx) const = 0;
    virtual matrix_t mkFunB(unsigned idx) const = 0;

    /** return or set vector of all factor matrices */
    virtual vector<matrix_t> mkFactorVec() const;
    virtual void mkFactorVec(vector<matrix_t> & v) const;

    virtual vector<matrix_t> mkFunAVec() const;
    virtual vector<matrix_t> mkFunBVec() const;

    /**     clear all submitted expectation counts */
    virtual void clearCounts() = 0;
    
    /** Return factor count */
    unsigned facCount() const {return facCount_;}

  protected:
    
    unsigned facCount_; // number of factors in set
  };


  class CompositeFactorSet : public AbstractBaseFactorSet {
  public:

    CompositeFactorSet(vector<AbsBasFacPtr_t> const factorPtrs) : AbstractBaseFactorSet( factorPtrs.size() ), factorPtrs_(factorPtrs) {}
    virtual ~CompositeFactorSet() {};

    /** submit expectations for matrix entries of each factor. */
    virtual void submitCounts(matrix_t const & counts, unsigned idx) {factorPtrs_[idx]->submitCounts(counts);}
    using AbstractBaseFactorSet::submitCounts; // bringing other definitions into this name space (hidden otherwise)

    /** return value may optionally report on success or other aspects of optimization. */
    virtual int optimizeParameters() {int success = 1; for (unsigned i = 0; i < facCount_; i++) success *= factorPtrs_[i]->optimizeParameters(); return success;}

    /** return or set matrix defining factor idx  */
    virtual matrix_t mkFactor(unsigned idx) const {return factorPtrs_[idx]->mkFactor();}
    virtual void mkFactor(matrix_t & m, unsigned idx) const  {return factorPtrs_[idx]->mkFactor(m);}
    virtual matrix_t mkFunA(unsigned idx) const {return factorPtrs_[idx]->mkFunA();}
    virtual matrix_t mkFunB(unsigned idx) const {return factorPtrs_[idx]->mkFunB();}

    using AbstractBaseFactorSet::mkFactor; // bringing other mkFactor definition into this name space (hidden otherwise)

    /** clear all submitted expectation counts */
    virtual void clearCounts() {for (unsigned i = 0; i < facCount_; i++) factorPtrs_[i]->clearCounts();}

  protected:

    vector<AbsBasFacPtr_t> const factorPtrs_; // smart pointers, defined above
  };


} // end namespace phy

#endif  // __Factors_h
