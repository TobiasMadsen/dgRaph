/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "Observations.h"
#include <sstream>

namespace phy {


  /****************************************
   TODO Split into another header
   ****************************************/

  /** Implementations of StateMaps following the PIMPL idiom **/
  class StateMapImplContinuous : public StateMapImpl {
  public:
    StateMapImplContinuous( unsigned bins, number_t minv, number_t maxv) : bins_(bins), minv_(minv), maxv_(maxv), convertIndexNumber_(minv,maxv,bins) { }
    
    virtual state_t const symbol2State(symbol_t const & s) const;
    virtual symbol_t state2Symbol(state_t i) const;
    virtual unsigned stateCount() const { return bins_; }
    virtual unsigned metaStateCount() const { return bins_; }
    virtual state_t symbolSize() const{ return 1;}
    virtual vector<symbol_t> const degeneracyVector(symbol_t const & s) const;
    virtual void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const;
    virtual void serialize(ostream & str) const;
  private:
    unsigned bins_;
    number_t minv_;
    number_t maxv_;
    vector<stateMask_t> metaState2StateMask_; //TODO make dynamic
    ConvertIndexNumber convertIndexNumber_;
  };


  class StateMapImplSymbol : public StateMapImpl {
  public:
    StateMapImplSymbol(string const & symbols) : state2Symbol_( stringToVectorOfStrings(symbols) ) { init(); }
    StateMapImplSymbol(vector<symbol_t> const & symbols) : state2Symbol_(symbols) { init(); }
    StateMapImplSymbol(vector<symbol_t> const & symbols, boost::unordered_map<symbol_t, vector<symbol_t> > const & metaSymbolDegeneracyMap) : state2Symbol_(symbols), degeneracyMap_(metaSymbolDegeneracyMap){ init(); }
    StateMapImplSymbol(StateMap const & staMap, unsigned n);

    virtual state_t const symbol2State(symbol_t const & s) const;
    virtual symbol_t state2Symbol(state_t i) const { return state2Symbol_[i]; }
    virtual unsigned stateCount() const { return stateCount_; } //Move to StateMapImpl?
    virtual unsigned metaStateCount() const { return metaStateCount_; }
    virtual state_t symbolSize() const { return symbolSize_; } //TODO Move to StateMapImpl
    virtual vector<symbol_t> const degeneracyVector(symbol_t const & s) const;
    virtual void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask_) const;
    virtual void serialize(ostream & str) const;
  private:
    void init();

    vector<symbol_t> state2Symbol_;
    boost::unordered_map<symbol_t, state_t> symbol2State_;
    boost::unordered_map<symbol_t, vector<symbol_t> > degeneracyMap_;
    unsigned stateCount_;
    unsigned metaStateCount_;
    unsigned symbolSize_;
  };

  class StateMapImplCount : public StateMapImpl {
  public:
    /** Constructor */
    StateMapImplCount(unsigned minv, unsigned maxv) : minv_(minv), maxv_(maxv) { }

    virtual state_t const symbol2State( symbol_t const & s) const;
    virtual symbol_t state2Symbol( state_t i) const;
    virtual unsigned stateCount() const { return maxv_-minv_+1;}
    virtual unsigned metaStateCount() const { return maxv_-minv_+1;}
    virtual state_t symbolSize() const{ return 1;}
    virtual vector<symbol_t> const degeneracyVector(symbol_t const & s) const;
    virtual void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const;
    virtual void serialize(ostream & str) const;

  private:
    unsigned minv_;
    unsigned maxv_;
    vector<stateMask_t> metaState2StateMask_; //TODO make dynamic
  };

  StateMapImplSymbol::StateMapImplSymbol(StateMap const & staMap, unsigned n){
    state2Symbol_ = mkMultiStateSymbols(staMap, n);
    boost::unordered_map<symbol_t, vector<symbol_t> > degMap;
    for (unsigned i = staMap.stateCount(); i < staMap.metaStateCount(); i++) {
      symbol_t sym = staMap.state2Symbol(i);
      degMap[sym] = staMap.degeneracyVector(sym);
    }
    degeneracyMap_ = mkMultiStateSymbolDegeneracyMap(degMap, staMap, n);

    init();
  }

  void StateMapImplSymbol::init(){
    stateCount_ = state2Symbol_.size();
    // add basic symbols to degeneracy map
    BOOST_FOREACH(symbol_t const & sym,  state2Symbol_)
      degeneracyMap_[sym] = vector<symbol_t>(1, sym);
    // add (only) metaSymbols to state2Symbols_
    symbol_t sym;
    vector<symbol_t> degSymVec;
    BOOST_FOREACH(boost::tie(sym, degSymVec), degeneracyMap_)
      if ( find(state2Symbol_.begin(), state2Symbol_.end(), sym) == state2Symbol_.end() ) // not found
	state2Symbol_.push_back(sym);
    metaStateCount_ = state2Symbol_.size();

    // symbol2State
    for (unsigned i = 0; i < metaStateCount_; i++)
      symbol2State_[ state2Symbol_[i] ] = i;
    //symbolSize
    symbolSize_= (stateCount_ > 0) ? state2Symbol_[0].size() : 0;
  }

  void StateMapImplSymbol::setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const{
    //Checks
    if(metaState2StateMask.size() != metaStateCount_)
      errorAbort("StateMapImplSymbol: metaState2StateMask_.size() != metaStateCount_");

    for (unsigned i = 0; i < metaStateCount_; i++) {
      symbol_t sym = state2Symbol(i);
      vector<symbol_t> const & v = degeneracyVector(sym);
      for (vector<symbol_t>::const_iterator it = v.begin(); it < v.end(); it++) {
	state_t j = symbol2State(*it);
	if (j >= stateCount_)
	  errorAbort("StateMaskMap::StateMaskMap: Degenerate symbol '" + sym + "' is defined in terms of another meta-symbol '" + *it + "'.");
	metaState2StateMask.at(i)(j) = true;
      }
    }
  }

  void StateMapImplSymbol::serialize(ostream & str) const {
    // symbols
    str << "SYMBOLS:\t";
    for (unsigned i = 0; i < stateCount_; i++)
      str << " " << state2Symbol(i);
    str << endl;

    // meta symbols
    str << "META_SYMBOLS:\t";
    for (unsigned i = stateCount_; i < metaStateCount_; i++) {
      vector<symbol_t> degVec = degeneracyVector( state2Symbol(i) );
      str << " " << state2Symbol(i) << " =";
      BOOST_FOREACH(string const & s, degVec)
	str << " " << s;
      str << ";";
    }
    str << endl;
  }

  void StateMapImplContinuous::setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const {
    if(metaState2StateMask.size() != stateCount())
      errorAbort("StateMapImplSymbol: metaState2StateMask_.size() != stateCount_");

    for(int i = 0 ; i < stateCount() ; ++i ){
      metaState2StateMask.at(i)(i) = true;
    }    
  }

  state_t const StateMapImplSymbol::symbol2State(symbol_t const & s) const {
    boost::unordered_map<symbol_t, state_t >::const_iterator it = symbol2State_.find(s);
    if ( it == symbol2State_.end() )
      errorAbort("StateMap::symbol2State: Symbol '" + s + "' not found in stateMap.");
    return it->second;
  }

  state_t const StateMapImplContinuous::symbol2State(symbol_t const & s) const {
    double sym;
    try{
      sym = boost::lexical_cast<double>(s);
    }
    catch( boost::bad_lexical_cast &){
      if( s == "NA" ){
	errorAbort("TODO: Implement NA for continuous observations");
      }
      errorAbort("Statemap::symbol2State: Symbol '" + s + "' could not be converted to double");
    }
    if(sym < minv_ or sym > maxv_) 
      errorAbort("StateMapImplContinuous: Symbol," + s + ", out of range");

    return bins_*(sym-minv_)/(maxv_-minv_);
  }

  symbol_t StateMapImplContinuous::state2Symbol(state_t i) const {
    std::stringstream s;
    s << convertIndexNumber_.indexToNumber(i);
    return s.str();
  }

  vector<symbol_t> const StateMapImplContinuous::degeneracyVector(symbol_t const & s) const{
    //TODO Refactor this function away
    errorAbort("This function is currently unavailable for continuous statemaps. Did you try to generate a multisymbol map from a continuous statemap?");
  }

  vector<symbol_t> const StateMapImplSymbol::degeneracyVector(symbol_t const & s) const{
    boost::unordered_map<symbol_t, vector<symbol_t> >::const_iterator it = degeneracyMap_.find(s);
    if ( it == degeneracyMap_.end() )
      errorAbort("StateMap::degeneracyVector: Symbol '" + s + "' not found in stateMap.");
    return it->second;
  }

  void StateMapImplContinuous::serialize(ostream & str) const{
    str << "REAL:\tTRUE" << std::endl;
    str << "BINS:\t" << bins_ << std::endl;
    str << "MIN:\t" << minv_ << std::endl;
    str << "MAX:\t" << maxv_ << std::endl;
  }

  state_t const StateMapImplCount::symbol2State( symbol_t const & s) const{
    int sym;
    try{
      sym = boost::lexical_cast<int>(s);
    }
    catch( boost::bad_lexical_cast &){
      if( s == "NA")
	errorAbort("TODO: Implement NA for count observations");
      errorAbort("StateMapImplCount::symbol2State: Symbol '" + s + "' could not be converted to integer");
    }
    if(sym < minv_ or sym > maxv_)
      errorAbort("StateMapImplCount::symbol2State: Symbol '" + s + "' out of range");
    return sym-minv_;
  }

  symbol_t StateMapImplCount::state2Symbol( state_t i) const{
    std::stringstream s;
    s << (minv_ + i);
    return s.str();
  }

  vector<symbol_t> const StateMapImplCount::degeneracyVector(symbol_t const & s) const{
    //TODO Refactor this function away!
    errorAbort("This function is currently unavailable or continuous statemaps. Did you try to generate a multisymbol map from a count statemap?");
  }

  void StateMapImplCount::setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const {
    if(metaState2StateMask.size() != stateCount())
      errorAbort("StateMapImplCount::setMetaState2StateMask: metaState2StateMask_.size() != stateCount_");

    for(int i = 0 ; i < stateCount() ; ++i )
      metaState2StateMask.at(i)(i) = true;
  }

  void StateMapImplCount::serialize(ostream & str) const {
    str << "REAL:\tCOUNT" << std::endl;
    str << "MIN:\t" << minv_ << std::endl;
    str << "MAX:\t" << maxv_ << std::endl;
  }


  /*****************************************
     StateMap stuff
   *****************************************/

  StateMap const & StateMap::operator=(StateMap const &rhs)
  {
    name_ = rhs.name_;
    pImpl = rhs.pImpl;
    return *this;
  }

  /* StateMap constructors */
  /** Constructor */
  StateMap::StateMap(string const & symbols, string const & name) : name_(name), pImpl(new StateMapImplSymbol(symbols)) { }

  /** Constructor */
  StateMap::StateMap(vector<symbol_t> const & symbols, string const & name) : name_(name), pImpl(new StateMapImplSymbol(symbols)) { }

  /** Constructor */
  StateMap::StateMap(vector<symbol_t> const & symbols, boost::unordered_map<symbol_t, vector<symbol_t> > const & metaSymbolDegeneracyMap, string const & name) : name_(name), pImpl(new StateMapImplSymbol(symbols, metaSymbolDegeneracyMap)) { }
 
  /** Constructor for continuous type StateMap */
  //  StateMap::StateMap( vector_t const & breakpoints, string const & name) : name_(name), pImpl(new StateMapImplContinuous(breakpoints)) { }

  /** Constructor for continuous type StateMap */
  StateMap::StateMap(unsigned bins, number_t minv, number_t maxv, string const & name) : name_(name), pImpl(new StateMapImplContinuous(bins, minv, maxv)) { }

  /** Constructor for count type StateMap */
  StateMap::StateMap(unsigned minv, unsigned maxv, string const & name) : name_(name), pImpl(new StateMapImplCount(minv,maxv)) { }

  StateMap::StateMap(StateMap const & staMap, unsigned n, string const & explicitName) : name_(explicitName), pImpl( new StateMapImplSymbol(staMap, n) )
  {
    if (name_.size() == 0)
      if ( staMap.name().size() )
	name_ = toString(n) + "-" + staMap.name();
  }

  vector<symbol_t> const StateMap::state2Symbol(vector<state_t> v) const
  {
    vector<symbol_t> u;
    u.reserve( v.size() );
    for (unsigned i = 0; i < v.size(); i++) {
      u.push_back( state2Symbol(v[i]) );
    }
    return u;
  }

  vector<state_t> const StateMap::symbol2State(vector<symbol_t> v) const
  {
    vector<state_t> u( v.size() );
    for (unsigned i = 0; i < v.size(); i++) {
      u.push_back( symbol2State(v[i]) );
    }
    return u;
  }

  string mkMultiStateSymbol(unsigned state, StateMap const & sm, unsigned n)
  {
    symbol_t symbol;
    unsigned base = sm.stateCount();
    for (int j = n - 1; j >= 0; j--) {
      unsigned val = power(base, j);
      state_t subState = state / val;
      symbol += sm.state2Symbol(subState);
      state = state % val;
    }
    return symbol;
  }


  // returns a symbol string of all possible n-long symbols constructed from the StateMap symbols
  vector<symbol_t> mkMultiStateSymbols(StateMap const & sm, unsigned n)
  {
    unsigned stateCount = power(sm.stateCount(), n);
    vector<symbol_t> symbols;
    symbols.reserve(stateCount);
    for (unsigned i = 0; i < stateCount; i++)
      symbols.push_back( mkMultiStateSymbol(i, sm, n) );
    return symbols;
  }
    

  void addBasicSymbolsToDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > & degMap, vector<symbol_t> const & symbols)
  {
    BOOST_FOREACH(symbol_t const & sym,  symbols)
      degMap[sym] = vector<symbol_t>(1, sym);
  }


  void addBasicSymbolsToDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > & degMap, StateMap const & staMap)
  {
    vector<symbol_t> symbols;
    for (unsigned i = 0; i < staMap.stateCount(); i++) 
      symbols.push_back( staMap.state2Symbol(i) );
    addBasicSymbolsToDegeneracyMap(degMap, symbols);
  }


  // take a degeneracy map and the corresponding state map and output a degeneracy map for the corresponding multi state symbol set.
  boost::unordered_map<symbol_t, vector<symbol_t> > mkMultiStateSymbolDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > const & degMap, StateMap const & staMap, unsigned n)
  {
    boost::unordered_map<symbol_t, vector<symbol_t> > extDegMap = degMap;
    addBasicSymbolsToDegeneracyMap(extDegMap, staMap);
    boost::unordered_map<symbol_t, vector<symbol_t> > mulDegMap = extDegMap;

    if (power(static_cast<unsigned int>(extDegMap.size()), static_cast<unsigned int>(n)) > 100000) // 
      errorAbort("Silly number of degenerate symbols in multi-symbol degeneracy map. Reimplement 'mkMultiStateDegeneracyMa' and related functions functions");

    // some loop internal variables
    boost::unordered_map<symbol_t, vector<symbol_t> > tmpMulDegMap;
    symbol_t preSym, appSym;
    vector<symbol_t> preDegSymVec, appDegSymVec;
    for (unsigned i = 1; i < n; i ++) {
      tmpMulDegMap.clear();
      BOOST_FOREACH(boost::tie(preSym, preDegSymVec), mulDegMap) {
	BOOST_FOREACH(boost::tie(appSym, appDegSymVec), extDegMap) {
	  symbol_t sym = preSym + appSym;
	  vector<symbol_t> degSymVec;
	  BOOST_FOREACH(symbol_t const & preDegSym, preDegSymVec) 
	    BOOST_FOREACH(symbol_t const & appDegSym, appDegSymVec) 
	    degSymVec.push_back(preDegSym + appDegSym);
	  tmpMulDegMap[sym] = degSymVec;
	}
      }
      mulDegMap = tmpMulDegMap;
    }
    return mulDegMap;
  }

  StateMap mkNucStateMap(void)
  {
    vector<symbol_t> symbols = stringToVectorOfStrings("ACGT");
    return StateMap(symbols, "nuc");
  }


  StateMap mkMetaNucStateMap(void)
  {
    vector<symbol_t> symbols = stringToVectorOfStrings("ACGT");
    boost::unordered_map<symbol_t, vector<symbol_t> > degMap;
    
    degMap["U"] = stringToVectorOfStrings("T"); 
    degMap["R"] = stringToVectorOfStrings("AG");
    degMap["Y"] = stringToVectorOfStrings("CT");
    degMap["K"] = stringToVectorOfStrings("GT");   
    degMap["M"] = stringToVectorOfStrings("CA");   
    degMap["S"] = stringToVectorOfStrings("CG");   
    degMap["W"] = stringToVectorOfStrings("AT");   
    degMap["B"] = stringToVectorOfStrings("CGT");  
    degMap["D"] = stringToVectorOfStrings("AGT");  
    degMap["H"] = stringToVectorOfStrings("ACT");  
    degMap["V"] = stringToVectorOfStrings("ACG");  
    degMap["N"] = stringToVectorOfStrings("ACGT");
    degMap["-"] = stringToVectorOfStrings("ACGT");
    degMap["."] = stringToVectorOfStrings("ACGT");

    return StateMap(symbols, degMap, "metaNuc");
  }


  // creates stateMap where each symbol consists of n nucleotides.
  StateMap mkMultiMetaNucStateMap(unsigned n)
  {
    StateMap metaNucMap = mkMetaNucStateMap();
    return StateMap(metaNucMap, n);
  }


  // creates stateMap where each symbol consists of n nucleotides.
  StateMap mkMultiNucStateMap(unsigned n)
  {
    StateMap nucMap = mkNucStateMap();
    return StateMap(nucMap, n);
  }


  StateMaskMap::StateMaskMap(StateMap const & staMap)
    : metaState2StateMask_( staMap.metaStateCount() , stateMask_t(staMap.stateCount(), false) )      
  {
    staMap.setMetaState2StateMask(metaState2StateMask_);
  }


  StateMaskMapSet::StateMaskMapSet(StateMap const & staMap, unsigned n) : stateMapCount_(n), stateMapVector_(n), stateMaskMapVector_(n)
  {
    StateMapPtr_t staMapPtr(new StateMap(staMap));
    StateMaskMapPtr_t staMasPtr(new StateMaskMap(staMap));

    for (unsigned i = 0; i < n ;i++) {
      stateMapVector_[i] = staMapPtr;
      stateMaskMapVector_[i] = staMasPtr;
    }
  }


  StateMaskMapSet::StateMaskMapSet(vector<StateMapPtr_t> const & stateMapVector) : stateMapCount_( stateMapVector.size() ), stateMapVector_(stateMapVector), stateMaskMapVector_(stateMapVector.size() ) 
  {
    unsigned n = stateMapCount_;
    for (unsigned i = 0; i < n ;i++) 
      if (stateMapVector_[i] != NULL)
	stateMaskMapVector_[i] = StateMaskMapPtr_t(new StateMaskMap(*stateMapVector_[i]));
      else {
	stateMaskMapVector_[i] = StateMaskMapPtr_t();
	//errorAbort("from StateMaskMapSet::StateMaskMapSet: Null StateMapPtr for index " + toString(i) + ". This should not happen. ");
      }
  }

  void StateMaskMapSet::symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec) const
  {
    assert(symVec.size() == stateMapCount_);

    for (unsigned i = 0; i < stateMapCount_; i++) {
      unsigned state = stateMapVector_[i]->symbol2State( symVec[i] );
      obsMasVec[i] = & stateMaskMapVector_[i]->metaState2StateMask(state);
    }
  }


  stateMaskVec_t StateMaskMapSet::symbols2StateMasks(vector<symbol_t> const & symVec) const
  {
    stateMaskVec_t obsMasVec(stateMapCount_);
    symbols2StateMasks(obsMasVec, symVec);
    return obsMasVec;
  }


  void StateMaskMapSet::symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec, vector<unsigned> const & seqToVarMap) const
  {
    assert( symVec.size() == seqToVarMap.size() );
    assert( obsMasVec.size() == stateMapCount_ );

    for (unsigned i = 0; i < stateMapCount_; i++)
      obsMasVec[i] = NULL;
    for (unsigned i = 0; i < symVec.size(); i++) {
      unsigned const & ranVarIdx = seqToVarMap[i];
      state_t state = stateMapVector_[ranVarIdx]->symbol2State( symVec[i] );
      obsMasVec[ ranVarIdx ] = & stateMaskMapVector_[ranVarIdx]->metaState2StateMask(state);
    }
  }

  unsigned long symbolCount(vector<string> const & strVec, unsigned offSet)
  {
    assert(offSet > 0);
    long seqSize = (strVec.size() > 0) ? strVec[0].size() : 0;
    unsigned long symCount = (seqSize + (offSet - 1) ) / offSet;  // number of symbols in input sequences. The result should be ceil for integer division.
    return symCount;
  }


  long unsigned symbolCount(long seqSize, unsigned offSet)
  {
    assert(offSet > 0);
    long unsigned symCount = (seqSize + (offSet - 1) ) / offSet;  // number of symbols in input sequences. The result should be ceil for integer division.
    return symCount;
  }

} // end namespace phy
