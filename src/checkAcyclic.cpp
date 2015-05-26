#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "rToCpp.h"

using namespace Rcpp;

bool visitVariable(int current, int parent, std::vector<std::vector<unsigned> > & facNbs, std::vector<std::vector<unsigned> > & varNbs, std::vector<bool> & facVisited, std::vector<bool> & varVisited);

bool visitFactor(int current, int parent, std::vector<std::vector<unsigned> > & facNbs, std::vector<std::vector<unsigned> > & varNbs, std::vector<bool> & facVisited, std::vector<bool> & varVisited);

bool visitVariable(int current, int parent, std::vector<std::vector<unsigned> > & facNbs, std::vector<std::vector<unsigned> > & varNbs, std::vector<bool> & facVisited, std::vector<bool> & varVisited){
  // Check if already visited
  if( varVisited.at(current) )
    return false;
  varVisited.at(current) = true;
  
  // Visit neighbors except for parent
  for(int i = 0; i < varNbs.at(current).size(); ++i){
    unsigned toVisit = varNbs.at(current).at(i);
    if(toVisit != parent)
      if( ! visitFactor(toVisit, current, facNbs, varNbs, facVisited, varVisited) )
	return false;
  }

  return true;
}

bool visitFactor(int current, int parent, std::vector<std::vector<unsigned> > & facNbs, std::vector<std::vector<unsigned> > & varNbs, std::vector<bool> & facVisited, std::vector<bool> & varVisited){
  // check if already visited
  if( facVisited.at(current) )
    return false;
  facVisited.at(current) = true;

  // Visit neighbors except for parent
  for(int i = 0; i < facNbs.at(current).size(); ++i){
    unsigned toVisit = facNbs.at(current).at(i);
    if(toVisit != parent)
      if( ! visitVariable(toVisit, current, facNbs, varNbs, facVisited, varVisited))
	  return false;
  }
    
  return true;
};

// [[Rcpp::export]]
bool checkAcyclic(List rFacNbs) {
   // Convert to vector of vectors
   std::vector<std::vector<unsigned> > facNbs = rNbsToNbs( rFacNbs);
   
   // Create vector of varNbs
   std::vector<std::vector<unsigned> > varNbs(0);
   for(int i = 0; i < facNbs.size(); ++i){
     for(int j = 0; j < facNbs.at(i).size(); ++j){
       varNbs.resize( (facNbs.at(i).at(j)+1 > varNbs.size())? 
		      (facNbs.at(i).at(j)+1):
		      varNbs.size() );
       varNbs.at( facNbs.at(i).at(j) ).push_back(i);
     }
   }
   
   // DFS
   std::vector<bool> facVisited(facNbs.size());
   std::vector<bool> varVisited(varNbs.size());

   for(int i = 0; i < facNbs.size(); ++i){
     if( ! facVisited.at(i) )
       if( ! visitFactor(i, -1, facNbs, varNbs, facVisited, varVisited) )
	 return false;
   }

   return true;
}
