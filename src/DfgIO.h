/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __DfgIO_h
#define __DfgIO_h

#include "DiscreteFactorGraph.h"
#include "Factors.h"
#include "Observations.h"

namespace phy {

  using namespace std;

  /** Structure holding a discrete factor graph (dfg) and convenient information for using it. */
  class DfgInfo 
  {
  public:
    /* Constructors **/
    DfgInfo(vector<string> const & varNames, vector<string> const & facNames, vector<AbsBasFacPtr_t> const & facVec, vector<StateMapPtr_t> const & stateMapVec, DFG const & dfg) :
      varNames(varNames), facNames(facNames), facVec(facVec), facSet(facVec), stateMapVec(stateMapVec), stateMaskMapSet(stateMapVec), dfg(dfg) {};

    DfgInfo(vector<string> const & varNames, 
	    vector<string> const & facNames, 
	    vector<string> const & potNames, 
	    vector<vector<unsigned> > const & facNeighbors, 
	    map<string, StateMapPtr_t> const & smMap, 
	    map<string, AbsBasFacPtr_t> const & facMap, 
	    map<string, string> const & var2smMap);

    /* Data **/
    vector<string> varNames;           ///< defines enumeration of random variables
    vector<string> facNames;           ///< defines enumeration of factors
    vector<string> subNames;           ///< defines enumeration of subscribed variables
    vector<AbsBasFacPtr_t> facVec;     ///< factor potential (pointer) for each factor. Each factor also holds the potential name
    CompositeFactorSet facSet;         ///< data structure for holding and optimizing factor potentials
    vector<StateMapPtr_t> stateMapVec; ///< stateMap for each variable
    StateMaskMapSet stateMaskMapSet;   ///< defines mapping between symbols and stateMasks for all variables
    DFG dfg;                           ///< Discrete factor graph
    
    // For subscriptions
    vector<unsigned> subscriptionFacs;  ///< subset of facVec that subscribes to variables
    vector< vector<unsigned> > subscribedVars; ///< the variables that the subscriptionFactors subscribe to

    void writeInfo( ostream & str );

    void updateFactors(vector<symbol_t> const & varVec, vector<unsigned> const & varMap );
  };

  /** Helper functions declared here for debug purposes */
  vector<AbsBasFacPtr_t> mkFacVec(vector<string> const & potNames, map<string, AbsBasFacPtr_t> const & facMap);
  vector<StateMapPtr_t> mkStateMapVec(vector<string> const & varNames, map<string, string> const & var2smMap, map<string, StateMapPtr_t> const & smMap);
  vector<unsigned> mkVarDimensions(vector<StateMapPtr_t> const & smVec);
  map<string, StateMapPtr_t> smVecToSmMap(vector<StateMapPtr_t> const & smVec);
  vector<string> mkPotNames(vector<AbsBasFacPtr_t> const & facVec);
} // end namespace phy

#endif  // __DfgIO_h
