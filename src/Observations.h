/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __Observations_h
#define __Observations_h

#include "utils.h"
#include "PhyDef.h"
#include "boost/tuple/tuple.hpp"
#include "boost/foreach.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp> 
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>

using namespace std;

namespace phy {

  ////////////////////////////////////////////////////////////////
  // StateMaps: symbol_t -> state_t
  ////////////////////////////////////////////////////////////////


  /** Map input string (symbol_t) symbols to integer observation
      space. Degenerate observation symbols that map to several states
      are allowed. These symbols are mapped to metaStates with values
      greater than the basic states. Multiple symbols may map to the
      same set of states, however, they are assigned different
      metaStates internally.*/

  class StateMapImpl;

  class StateMapImpl{
  public:
    //    virtual ~StateMapImpl();
    virtual state_t const symbol2State(symbol_t const & s) const = 0;
    virtual symbol_t state2Symbol(state_t i) const = 0;
    //virtual symbol_t state2Symbol(state_t i, symbol_t & s) const = 0;


    /** returns degeneracy vector for symbol s */
    virtual vector<symbol_t> const degeneracyVector(symbol_t const & s) const = 0;

    /** returns number of basic states */
    virtual unsigned stateCount() const = 0;

    /** returns number of meta states (includes basic states) */
    virtual unsigned metaStateCount() const = 0;

    /** returns length (in chars) of symbols */
    virtual state_t symbolSize() const = 0;

    /** serialize statemap */
    virtual void serialize(ostream & os) const = 0;

    virtual void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask_) const = 0;
  };

  class StateMapImplContinuous;
  class StateMapImplSymbol;
  typedef boost::shared_ptr<StateMapImpl> StateMapImplPtr_t;
  
  class StateMap {
  public:
    
    /** Constructor */
    StateMap(string const & symbols, string const & name = "");

    /** Constructor */
    StateMap(vector<symbol_t> const & symbols, string const & name = "") ;

    /** Constructor */
    StateMap(vector<symbol_t> const & symbols, boost::unordered_map<symbol_t, vector<symbol_t> > const & metaSymbolDegeneracyMap, string const & name = "");
 
    /** Constructor of n-long multiSymbols based on a basic stateMap. If staMap has a name, and no explicit name is given, the name of the constructed StateMap will be n-name.*/
    StateMap(StateMap const & staMap, unsigned n, string const & explicitName = "");

    /** Constructor for continuous type StateMap*/
    StateMap(unsigned bins, number_t minv, number_t maxv, string const & name = "");

    /** Constructor for count type StateMap */
    StateMap(unsigned minv, unsigned maxv, string const & name = "");

    /** Copy assignment operator. */
    StateMap const & operator=(StateMap const &rhs);

    /** Returns symbol corresponding to state i */
    symbol_t state2Symbol(state_t i) const {return pImpl->state2Symbol(i);}
    vector<symbol_t> const state2Symbol(vector<state_t> v) const;

    /** Returns state corresponding to symbol s. Aborts on nonexisting symbols. */
    state_t const symbol2State(symbol_t const & s) const { return pImpl->symbol2State(s);}
    vector<state_t> const symbol2State(vector<symbol_t> v) const;

    /** returns degeneracy vector for symbol s */
    vector<symbol_t> const degeneracyVector(symbol_t const & s) const { return pImpl->degeneracyVector(s);}

    /** returns number of basic states */
    unsigned stateCount() const {return pImpl->stateCount();}

    /** returns number of meta states (includes basic states) */
    unsigned metaStateCount() const {return pImpl->metaStateCount();}

    /** returns length (in chars) of symbols */
    state_t symbolSize() const  {return pImpl->symbolSize(); }

    /** Set metaState2StateMask for the StateMaskMap class */
    void setMetaState2StateMask(vector<stateMask_t> & metaState2StateMask) const { pImpl->setMetaState2StateMask( metaState2StateMask);}

    /** Serialize statemap */
    void serialize(ostream &os) const {return pImpl->serialize( os);}


    /** Returns StateMap name. Name is empty ("") if not a canonical
	type. Useful in IO-functions. */
    string const & name() const {return name_;}

  protected:
    string name_;
  private:
    StateMapImplPtr_t pImpl;
  };

  /** returns an n-state symbol with state 'state' based on stateMap. */
  symbol_t mkMultiStateSymbol(unsigned state, StateMap const & sm, unsigned n);

  /** returns a vector of all n-state symbols composed from single state symbols. Used to create multi nucleotide stateMaps. */
  vector<symbol_t> mkMultiStateSymbols(StateMap const & sm, unsigned n);

  /** Construct the degeneracy map for n-state symbols baed on the degeneracy map given in the input StateMap (staMap). */
  boost::unordered_map<symbol_t, vector<symbol_t> > mkMultiStateSymbolDegeneracyMap(boost::unordered_map<symbol_t, vector<symbol_t> > const & degMap, StateMap const & staMap, unsigned n);


  ////////////////////////////////////////////////////////////////
  // Nucleotide specific StateMaps
  ////////////////////////////////////////////////////////////////
   
  /** Return a symbol map for basic nucleotides */
  StateMap mkNucStateMap();

  /** Return a symbol map for basic nucleotides and IUPAC symbols. Can be used with DNA as well as RNA. */
  StateMap mkMetaNucStateMap();

  /** Returns a stateMap where each symbol consists of n nucleotides. */
  StateMap mkMultiNucStateMap(unsigned n);

  /** Returns a stateMap where each symbol consists of n nucleotides. Include degenerate symbols and should not be used for more than n=3.*/
  StateMap mkMultiMetaNucStateMap(unsigned n);


  ////////////////////////////////////////////////////////////////
  // StateMaskMaps: state_t -> stateMask_t
  ////////////////////////////////////////////////////////////////

  /* StateMaskMap allows conversion from (meta)states to stateMask_t,
     which are vectors of bools denoting if a state is observed or
     not. A stateMask_t is the natural input data structure for, e.g.,
     discrete factor graphs. **/
  class StateMaskMap {
  public:
    //TODO: could allow dynamic calculation of stateMasks to save space Especially for continous factors!
    StateMaskMap(StateMap const & staMap);

    /** Return the bit vector corresponding to state i */
    stateMask_t const & metaState2StateMask(state_t i) const {return metaState2StateMask_.at(i);} 

  protected:
    //bool isCont_;
    vector<stateMask_t> metaState2StateMask_;
  };


  ////////////////////////////////////////////////////////////////
  // StateMaskMapSet: vector<symbol_t> -> stateMaskVec_t
  ////////////////////////////////////////////////////////////////

  typedef boost::shared_ptr<StateMaskMap> StateMaskMapPtr_t;
  typedef boost::shared_ptr<StateMap> StateMapPtr_t;

  /** Class that holds all the stateMasks for a given indexed set of
      random variables. Since the member functions return pointers to
      internally stored stateMasks, which are normally used in
      calculations on a factorGraph or similar, this datastructure
      must exist as long as calculations are performed or as long as
      the stateMasks are needed. */
  class StateMaskMapSet {
  public:

    /** Constructor to use when all random variables are over the same state space (i.e, all use the same StateMap). */
    StateMaskMapSet(StateMap const & staMap, unsigned n);

    /** Constructor to use when random variables are defined over different state spaces (i.e., they use different StateMaps). If StateMap is not defined for a random variable, a NULL pointer can be given.*/
    StateMaskMapSet(vector<StateMapPtr_t> const & stateMapVector);

    /** Convert a vector of symbols to a vector of observation
	masks. The indexing of symVec corresponds to the indexing of
	StateMaskMaps in StateMaskMapSet (which should correspond to the
	indexing of random variables in, e.g., a factor graph). The
	returned vector will have a stateMask pointer for each random
	variable. For unobserved random variables the version below
	can be used. */
    void symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec) const;
    stateMaskVec_t symbols2StateMasks(vector<symbol_t> const & symVec) const;

    /** Same as above but allows for unobserved variables. The randomVariableMap defines which random variable (and thus stateMask index) each index of symVec pertains to. obsMasVec[i] is set to null (missing data) if not pointed to by varMap. */
    stateMaskVec_t symbols2StateMasks(vector<symbol_t> const & symVec, vector<unsigned> const & varMap) const;
    void symbols2StateMasks(stateMaskVec_t & obsMasVec, vector<symbol_t> const & symVec, vector<unsigned> const & varMap) const;

    unsigned stateMapCount(void) const {return stateMapCount_;}
    unsigned symbolSize(unsigned idx) const {return stateMapVector_[idx]->symbolSize();}
    
  protected:
   
    unsigned stateMapCount_;
    vector<StateMapPtr_t> stateMapVector_;
    vector<StateMaskMapPtr_t> stateMaskMapVector_;
   };

  ////////////////////////////////////////////////////////////////
  // wrappers for making symbol_t and vector<symbol_t> from various
  // data types.
  // These are made inline for efficiency.
  ////////////////////////////////////////////////////////////////

  /** precondition: the symbol (sym) is of correct size!
      the symbol start position (symStartPos) is the first position of the symbol, I.e., largest meaningful pos is seqSize. */
  inline void mkSymbol(symbol_t & sym, unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar = '.');
  inline symbol_t mkSymbol(unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar = '.');
  
  inline void mkSymbolVector(vector<symbol_t> & symVec, vector<unsigned> const & symSizeVec, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData = '.');
  inline void mkSymbolVector(vector<symbol_t> & symVec, unsigned const commonSymSize, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData = '.');

  // A diSymbol is composed of left and a right part.
  // precondition: the symbol (sym) is of correct size!
  // the symbol start position (symStartPos) is the first position of the symbol, I.e., largest meaningful pos is seqSize.
  inline void mkDiSymbol(symbol_t & sym, unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar = '.');
  inline symbol_t mkDiSymbol(unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar = '.');

  // A diSymbol is composed of a left and a right part.
  inline void mkDiSymbolVector(vector<symbol_t> & symVec, unsigned const commonSizeLeft, unsigned const commonSizeRight, vector<string> const & strVec, unsigned const seqSize, unsigned const leftStartPos, unsigned const rightStartPos, char missingDataChar = '.');

  /** make missing data symbol. Precondition: sym is of correct size. */
  inline void mkMissingDataSymbol(string & sym, unsigned symSize, char const missingDataChar = '.');

  ////////////////////////////////////////////////////////////////
  // Helper functions
  ////////////////////////////////////////////////////////////////

  /** return the number of symbols in SeqData given a certain fixed offset between each. */
  unsigned long symbolCount(vector<string> const & strVec, unsigned offSet);
  unsigned long symbolCount(long seqSize, unsigned offSet);

  ////////////////////////////////////////////////////////////////
  // inline function definitions
  ////////////////////////////////////////////////////////////////

  inline void mkSymbol(symbol_t & sym, unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar)
  {
    long pos;
      for (unsigned i = 0; i < symSize; i++) {
	pos = symStartPos + i;
      if (pos < 0)
	sym[i] = missingDataChar;
      else if (pos >= seqSize)
	sym[i] = missingDataChar;
      else 
	sym[i] = seq[pos];
    }
  }


  inline symbol_t mkSymbol(unsigned symSize, string const & seq, long seqSize, long symStartPos, char missingDataChar)
  {
    symbol_t sym(symSize, missingDataChar);
    mkSymbol(sym, symSize, seq, seqSize, symStartPos, missingDataChar);
    return sym;
  }
  
  
  inline void mkSymbolVector(vector<symbol_t> & symVec, vector<unsigned> const & symSizeVec, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData)
  {
    for (unsigned i = 0; i < symVec.size(); i++)
      mkSymbol(symVec[i], symSizeVec[i], strVec[i], seqSize, symStartPos, missingData);
  }


  inline void mkSymbolVector(vector<symbol_t> & symVec, unsigned const commonSymSize, vector<string> const & strVec, unsigned seqSize, unsigned symStartPos, char missingData)
  {
    for (unsigned i = 0; i < symVec.size(); i++)
      mkSymbol(symVec[i], commonSymSize, strVec[i], seqSize, symStartPos, missingData);
  }


  // A diSymbol is composed of left and a right part.
  // precondition: the symbol (sym) is of correct size!
  // the symbol start position (symStartPos) is the first position of the symbol, I.e., largest meaningful pos is seqSize.
  inline void mkDiSymbol(symbol_t & sym, unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar)
  {
    // left part
    long pos;
    for (unsigned i = 0; i < sizeLeft; i++) {
      pos = leftStartPos + i;
      if (pos < 0)
	sym[i] = missingDataChar;
      else if (pos >= seqSize)
	sym[i] = missingDataChar;
      else 
	sym[i] = seq[pos];
    }

    // right part
    for (unsigned i = 0; i < sizeRight; i++) {
      pos = rightStartPos + i;
      if (pos < 0)
	sym[sizeLeft + i] = missingDataChar;
      else if (pos >= seqSize)
	sym[sizeLeft + i] = missingDataChar;
      else 
	sym[sizeLeft + i] = seq[pos];
    }
  }


  inline symbol_t mkDiSymbol(unsigned const sizeLeft, unsigned const sizeRight, string const & seq, long const seqSize, long const leftStartPos, long const rightStartPos, char const missingDataChar)
  {
    symbol_t sym(sizeLeft + sizeRight, missingDataChar);
    mkDiSymbol(sym, sizeLeft, sizeRight, seq, seqSize, leftStartPos, rightStartPos, missingDataChar);
    return sym;
  }


  // A diSymbol is composed of a left and a right part.
  inline void mkDiSymbolVector(vector<symbol_t> & symVec, unsigned const commonSizeLeft, unsigned const commonSizeRight, vector<string> const & strVec, unsigned const seqSize, unsigned const leftStartPos, unsigned const rightStartPos, char missingDataChar)
  {
    for (unsigned i = 0; i < symVec.size(); i++)
      mkDiSymbol(symVec[i], commonSizeLeft, commonSizeRight, strVec[i], seqSize, leftStartPos, rightStartPos, missingDataChar);
  }

  inline void mkMissingDataSymbol(string & sym, unsigned symSize, char const missingDataChar)
  {
    for (unsigned i = 0; i < symSize; i++)
      sym[i] = missingDataChar;
  }


} // end namespace phy


#endif  //__Observations_h
