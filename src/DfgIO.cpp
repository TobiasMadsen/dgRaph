/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "DfgIO.h"
#include <Rcpp.h>

namespace phy {


  // helper functions, not declared in header
  vector<AbsBasFacPtr_t> mkFacVec(vector<string> const & potNames, map<string, AbsBasFacPtr_t> const & facMap)
  {
    vector<AbsBasFacPtr_t> facVec;
    BOOST_FOREACH(string const & pot, potNames) {
      map<string, AbsBasFacPtr_t>::const_iterator iter = facMap.find(pot);
      if ( iter == facMap.end() )
	errorAbort("From mkFacVec: factor with name '" + pot + "' not found in factor map.");
      facVec.push_back(iter->second);
    }
    return facVec;
  }


  vector<StateMapPtr_t> mkStateMapVec(vector<string> const & varNames, map<string, string> const & var2smMap, map<string, StateMapPtr_t> const & smMap)
  {
    vector<StateMapPtr_t> smVec;
    BOOST_FOREACH(string const & var, varNames) {
      map<string, string>::const_iterator iter1 = var2smMap.find(var);
      if ( iter1 == var2smMap.end() )
	errorAbort("From mkStateMapVec: variable with name '" + var + "' not found in variable to stateMap map.");
      string const sm = iter1->second;

      map<string, StateMapPtr_t>::const_iterator iter2 = smMap.find(sm);
      if ( iter2 == smMap.end() )
	errorAbort("From mkStateMapVec: state map with name '" + sm + "' not found in stateMap map.");
      smVec.push_back(iter2->second);
    }
    return smVec;
  }


  vector<unsigned> mkVarDimensions(vector<StateMapPtr_t> const & smVec)
  {
    vector<unsigned> varDim;
    BOOST_FOREACH(StateMapPtr_t const & smPtr, smVec)
      varDim.push_back( smPtr->stateCount() );
    return varDim;
  }


  map<string, StateMapPtr_t> smVecToSmMap(vector<StateMapPtr_t> const & smVec)
  {
    map<string, StateMapPtr_t> smMap;
    BOOST_FOREACH(StateMapPtr_t const & smPtr, smVec) {
      if (smPtr->name().size() == 0)
	errorAbort("from smVecToSmMap: attempt to add stateMap lacking name a to stateMap map");
      smMap[smPtr->name()] = smPtr;
    }
    return smMap;
  }


  vector<string> mkPotNames(vector<AbsBasFacPtr_t> const & facVec)
  {
    vector<string> potVec;
    BOOST_FOREACH(AbsBasFacPtr_t const & facPtr, facVec)
      potVec.push_back( facPtr->name() );
    return potVec;
  }

  // definition of header declared functions
  DfgInfo::DfgInfo(vector<string> const & varNames, 
		   vector<string> const & facNames, 
		   vector<string> const & potNames, 
		   vector<vector<unsigned> > const & facNeighbors, 
		   map<string, StateMapPtr_t> const & smMap, 
		   map<string, AbsBasFacPtr_t> const & facMap, 
		   map<string, string> const & var2smMap) :
    varNames(varNames), 
    facNames(facNames), 
    facVec( mkFacVec(potNames, facMap) ),
    facSet(facVec),
    stateMapVec( mkStateMapVec(varNames, var2smMap, smMap) ), 
    stateMaskMapSet(stateMapVec),
    dfg( mkVarDimensions(stateMapVec), facSet.mkFactorVec(), facNeighbors )
  {
    //Add subscribed variables to subNames
    for(int facIdx = 0; facIdx < facVec.size(); ++facIdx){
      vector<string> sn = facVec.at(facIdx)->getSubscriptions();
      if( sn.size() == 0 )
	continue;

      // List of factors that subscribe to at least one variable
      subscriptionFacs.push_back(facIdx);
      subscribedVars.push_back( vector<unsigned>());
      
      //For each element in sn check if already present otherwise insert
      for(int subIdx = 0; subIdx < sn.size(); ++subIdx){
	vector<string>::iterator it = std::find( subNames.begin(), subNames.end(), sn.at(subIdx));
	subscribedVars.back().push_back( it - subNames.begin() );
	if( it == subNames.end() )
	  subNames.push_back( sn.at(subIdx));
      }

    }

  };

  void DfgInfo::writeInfo(ostream & str )
  {
    for(unsigned i = 0; i < facVec.size(); ++i){
      str << "NAME:\t" << facVec[i]->name() << endl;
      facVec[i]->serialize(str);
    }
  }

  // subNames enumerates subscription variables
  // varVec contains symbols for each subscription variable
  // varMap[i] gives the index in varVec corresponding to subNames[i] (can be obtained using invMap in NamedData)
  void DfgInfo::updateFactors(vector<symbol_t> const & varVec, vector<unsigned> const & varMap){
    vector<matrix_t> facPot( subscriptionFacs.size() );
    for(int f = 0; f < subscriptionFacs.size(); ++f){
      //Call update with correct variables
      vector<string> vars;

      //Fill the correct variables into vars 
      for(int i = 0; i < subscribedVars.at(f).size(); ++i){
	unsigned varPos;
	symbol_t var;
	try{
	  varPos = varMap.at(subscribedVars.at(f).at(i));
	}
	catch(const std::out_of_range& e){
	  errorAbort("Out of range  DfgIO updateFactors 1");
	}
	try{
	  var = varVec.at(varPos);
	}
	catch(const std::out_of_range& e){
	  Rcpp::Rcout << "Following varpos was not found in vector: " << varMap.at(subscribedVars.at(f).at(i)) << std::endl;
	  errorAbort("Out of range  DfgIO updateFactors 2");
	}
	vars.push_back( var); 			  
      }
      facVec.at(subscriptionFacs.at(f))->update(vars);
      facPot.at(f) = facVec.at(subscriptionFacs.at(f))->mkFactor();
    }
    dfg.resetFactorPotentials( facPot, subscriptionFacs );
  }

  map<string, vector<string> > mkStateMap2varNamesMap(vector<string> const & varNames, vector<StateMapPtr_t> const & stateMapVec)
  {
    assert( varNames.size() == stateMapVec.size() );
    map<string, vector<string> > sm2varNames;
    for (unsigned i = 0; i < varNames.size(); i++) {
      string smName  = stateMapVec[i]->name();
      string varName = varNames[i];
      map<string, vector<string> >::const_iterator iter=sm2varNames.find(smName);
      if ( iter == sm2varNames.end() ) // not found
	sm2varNames[smName] = vector<string>();
      sm2varNames[smName].push_back(varName); 
    }
    return sm2varNames;
  }

  vector<unsigned> mkFacNbs(vector<string> & varNames, vector<string> facNbsNames)
  // facNbsNames is an (ordered) list of random variables neighboring a specific factor
  // new variable names will be added to varNames and the facNbsNames converted into a list of varNames indices
  {
    vector<unsigned> facNbs;
    // add var names if not present
    BOOST_FOREACH(string const & var, facNbsNames) {
      vector<string>::const_iterator iter = find(varNames.begin(), varNames.end(), var);
      if ( iter == varNames.end() )
	varNames.push_back(var);
      facNbs.push_back( getIndex(varNames, var) );
    }
    return facNbs;
  }

} // end namespace phy
