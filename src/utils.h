/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __utils_h
#define __utils_h

////////////////////////////////////////////////////////////////
// This header defines convenient text and string manipulation
// functions
//
////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>
#include <fstream>
#include <utility>
#include <algorithm>
#include "Definitions.h"

#ifndef BOOST_TEST
#include <Rcpp.h>
#endif

namespace dgRaph {

  // ERROR CALLS

  /** Output message and terminate. */
  inline void errorAbort(string const & message = "");


  // CONVERSION 

  /** Returns a string representation vector v where each element is separated by sep. Use sep="" for concatenation. */
  template <class T>
  inline string toSepString (vector<T> const & v, string const & sep = ",", unsigned prec = 0);

  /** Returns the string representation of t. */
  template <class T>
  inline string toString(T const & t, unsigned prec = 0);

  /** Returns a string representation vector v. */
  template <class T>
  inline string toString (vector<T> const & v, unsigned prec = 0);

  /** Returns the string representation of pair p. */
  template <class T, class U>
  inline string toString (pair<T, U> const & p, unsigned prec = 0);

  /** Convert string s to val. */
  template <class T>
  inline void fromString(string const & s, T & val);

  /** Convert vector of strings v to vector of values u. */
  template <class T>
  inline void fromString(vector<string> const & v, vector<T> & u);



  // STRING MANIPULATIONS

  string strip(string const & s, string const & sep = " \t\n");

  vector<string> split(string const & s, string const & sep = " \t\n");
  inline vector<string> stringToVectorOfStrings(string const & s);

  // MATH RELATED

  /** Integer power function. Lacking from math.h */
  inline unsigned power(unsigned x, unsigned p);
  /** For consistency with NTL */
  inline double power(double x, double p) {return pow(x, p);}
  inline double power(double x, unsigned p) {return pow(x, static_cast<int>(p) );}

  class ConvertIndexNumber {
  public:
    ConvertIndexNumber(number_t const & minv, number_t const & maxv, unsigned const & bins) : minv_(minv), maxv_(maxv), bins_(bins) { 
      toNumberAlpha_ = (maxv_-minv_)/bins_;
      toNumberBeta_ = minv_+(maxv_-minv_)/2/bins_;
      toIndexAlpha_ = (maxv_-minv_)/bins_;
      toIndexBeta_ = -minv_*(maxv_-minv_)/bins_;
    }

    number_t const & getToNumberAlpha() const{ return toNumberAlpha_;}
    number_t const & getToNumberBeta() const{ return toNumberBeta_;}
    number_t const & getToIndexAlpha() const{ return toIndexAlpha_;}
    number_t const & getToIndexBeta() const{ return toIndexBeta_;}
    unsigned numberToIndex(number_t s) const {
      if( s < minv_ or s > maxv_ ) errorAbort("ConvertIndexNumber: NumberToIndex: Number out of range: Might change implementation return 0 if s < min and bins-1 if s > max");
      return bins_*(s-minv_)/(maxv_-minv_);
    }
    number_t indexToNumber(unsigned i) const {
      return i*toNumberAlpha_+toNumberBeta_;
    }

  private:
    number_t minv_;
    number_t maxv_;
    unsigned bins_;
    number_t toNumberAlpha_;
    number_t toNumberBeta_;
    number_t toIndexAlpha_;
    number_t toIndexBeta_;
  };

  /** Return accumulated log of values in v. */
  template<class V>
  number_t accLog(V const & v);


  // COMPARISON AND LOOK-UP

  /** Return the number at positions at which s1 and s2 differ. Precondition: s1 and s2 must be of the same length. */
  inline unsigned hammingDistance(string const & s1, string const & s2);
 
  /** Return true if element e is found in vector v. */
  template <class T>
  inline bool hasElement(vector<T> const & v, T const & e);

  template <class T>
  inline unsigned getIndex(vector<T> const & v, T const & e);

  /** Returns index map from subsetNames to allNames. */
  template <class T>
  vector<unsigned> mkSubsetMap(vector<T> const & allNames, vector<T> const & subsetNames);

  /** Return subset of v. */
  template <class T>
  void mkSubset(vector<T> const & v, vector<unsigned> const & subSetMap, vector<T> & subSet);  

  /** Same as above but returns result. */
  template <class T>
  vector<T> mkSubset(vector<T> const & v, vector<unsigned> const & subSetMap);


  // FILE PARSING

  /** parse parameters values from input strings of the type "tag val". The default sepation between tag and val is whitespace.*/
  template <class T>
  void getFeature(string inStr, string const & tag, T & val, string const & sep = " \t\n");

  /** Same as above, but parses a stream rather than a string.*/
  template <class T>
  void getFeature(istream & str, string const & tag, T & val);

  /** Open fileName and assign stream to f. Reports errors */
  inline void openInFile(ifstream & f, string const & fileName);

  /** Open fileName and assign stream to f. Reports errors */
  inline void openOutFile(ofstream & f, string const & fileName);

  /** Ignores characters in stream until (and including) the next newline */
  istream & skipLine(istream & str, unsigned maxChar = 1000000);

  /** Ignore all stream characters which are whitespace and the rest of all lines containing a '#' */
  istream & skipWhiteSpaceAndComments(istream & str);

  /** Return false if first character on line is whitespace. Skips comment lines (starting with '#') */
  bool moreTags(istream & str);

  /** Return the path prefix of filePath. Returns the empty string if no path prefix found. */
  string pathPrefix(string const & filePath);


  // INLINE IMPLEMTATIONS -- not part of interface

  inline void errorAbort(string const & message)
  {
    #ifdef BOOST_TEST
    std::cout << "Error. Program terminating. " << message << endl;
    exit(1);
    #else
    Rcpp::Rcout << "Error. Program terminating. " << message << endl;
    Rcpp::stop("Error");
    #endif
  }

  template <class T>
  inline string toString(T const & t, unsigned prec)
  {
    stringstream ss;
    if (prec)
      ss.precision(prec);
    ss << t;
    return ss.str();
  }


  // Returns a string representation vector v where each element is separated by sep. Use sep="" for concatenation.
  template <class T>
  inline string toSepString (vector<T> const & v, string const & sep, unsigned prec)
  {
    stringstream ss;
    if (prec)
      ss.precision(prec);
    unsigned size = v.size();
    for (unsigned i = 0; i < size; i++) {
      ss << v[i];
      if (i < size - 1)
	ss << sep;
    }
    return ss.str();
  }

  // Returns a string representation vector v. 
  template <class T>
  inline string toString (vector<T> const & v, unsigned prec)
  {
    stringstream ss;
    unsigned size = v.size();
    ss << "[" << size << "]" << "(" << toSepString(v, ",", prec) << ")";
    return ss.str();
  }


  // Returns a string representation pair p. 
  template <class T, class U>
  inline string toString (pair<T, U> const & p, unsigned prec)
  {
    stringstream ss;
    if (prec)
      ss.precision(prec);
    ss << "(" << p.first << ", " << p.second << ")";
    return ss.str();
  }


  // Convert string s to val. */
  template <class T>
  inline void fromString(string const & s, T & val)
  {
    istringstream ss;
    ss.str(s);
    ss >> val;
  }


  // Convert vector of strings v to vector of values u. 
  template <class T>
  inline void fromString(vector<string> const & v, vector<T> & u)
  {
    assert( v.size() == u.size() );
    for (unsigned i = 0; i < v.size(); i++)
      fromString(v[i], u[i]);
  }

  
  // Integer power function. Lacking from math.h */
  inline unsigned power(unsigned x, unsigned p) 
  {
    int val = 1;
    while(p--)
      val *= x;
    return val;
  }

  template<class V>
  number_t accLog(V const & v)
  {
    number_t r = 0;
    for (unsigned i = 0; i < v.size(); i ++)
      r += log(v[i]);
    return r;
  }


  inline vector<string> stringToVectorOfStrings(string const & s)
  {
    vector<string> v;
    v.reserve(s.size() );
    for (unsigned i = 0; i < s.size(); i ++)
      v.push_back(string(1, s[i]) );
    return v;
  }


  // Return the number at positions at which s1 and s2 differ. Precondition: s1 and s2 must be of the same length. */
  inline unsigned hammingDistance(string const & s1, string const & s2)
  {
    assert( s1.size() == s2.size() );
    unsigned n = 0;
    for (unsigned i = 0; i < s1.size(); i++)
      if (s1[i] != s2[i])
	n++;
    return n;
  }
  

  // Return true if element e is found in vector v. */
  template <class T>
  inline bool hasElement(vector<T> const & v, T const & e)
  {
    typename vector<T>::const_iterator result = find(v.begin(), v.end(), e);
    return ( result != v.end() );
  }
  

  template <class T>
  inline unsigned getIndex(vector<T> const & v, T const & e)
  {
    typename vector<T>::const_iterator result = find(v.begin(), v.end(), e);
    if ( result == v.end() )
      errorAbort("Element (" + toString(e) + ") not found in vector (printed below). \n" + toString(v) + "\n");
    
    return result - v.begin();
  }


  template <class T>
  vector<unsigned> mkSubsetMap(vector<T> const & allElements, vector<T> const & subsetElements)
  {
    vector<unsigned> subsetMap( subsetElements.size() );
    for (unsigned i = 0; i < subsetElements.size(); i++) {
      T const & element = subsetElements[i];
      typename vector<T>::const_iterator result = find(allElements.begin(), allElements.end(), element);
      if (result == allElements.end() )
	errorAbort("From mkSubsetMap: No element with value '" + string(element) + "' found.");
      subsetMap[i] = result - allElements.begin();
    }
    return subsetMap;
  }

  template <class T>
  vector<unsigned> mkMap(vector<T> const & toSet, vector<T> const & fromSet)
  {
    vector<unsigned> map( fromSet.size() );
    for (unsigned i = 0; i < fromSet.size(); i++) {
      T const & element = fromSet[i];
      typename vector<T>::const_iterator result = find(toSet.begin(), toSet.end(), element);
      if (result == toSet.end() )
	map[i] = -1;
      else
	map[i] = std::distance(toSet.begin(),result);
    }
    return map;
  }


  template <class T>
  void mkSubset(vector<T> const & v, vector<unsigned> const & subSetMap, vector<T> & subSet)
  {
    assert( subSetMap.size() == subSet.size() );
    for (unsigned i = 0; i < subSetMap.size(); i++) {
      unsigned idx = subSetMap[i];
      assert(idx < v.size());
      subSet[i] = v[idx];
    }
  }


  template <class T>
  vector<T> mkSubset(vector<T> const & v, vector<unsigned> const & subSetMap)
  {
    vector<T> subSet( subSetMap.size() );
    mkSubset(v, subSetMap, subSet);
    return subSet;
  }


  // parse parameters values from input strings of the type "tag val". The default sepation between tag and val is whitespace.*/
  template <class T>
  void getFeature(string inStr, string const & tag, T & val, string const & sep)
  {
    string::size_type pos = inStr.find(tag);
    if (pos == string::npos)
      errorAbort("The tag:\n  '" + tag + "'\nwas not found in input string:\n  '" + inStr + "'.");
    inStr.erase(0, pos + tag.size() );
    fromString(inStr, val);
  }

  // Same as above, but parses a stream rather than a string.
  // Will skip rest of line
  template <class T>
  void getFeature(istream & str, string const & tag, T & val)
  {
    string s;
    str >> s;
    if (s != tag)
      errorAbort("The tag found '" + s + "' does not match the expected tag '" + tag + "'.");
    str >> val;
  }

  // Same as above, will skip rest of line
  template <class T>
  void getFeatureAndSkipLine(istream & str, string const & tag, T & val)
  {
    getFeature(str, tag, val);
    skipLine(str);
  }

  // Get a tag and corresponding value as string
  // To be used if multiple tags can be used for instance when we have different parameterizations
  inline void getTagFeatureAndSkipLine(istream & str, string & tag, string & val)
  {
    str >> tag;
    str >> val;
    skipLine(str);
  }

  

  inline void openInFile(ifstream & f, string const & fileName)
  {
    f.open( fileName.c_str() );
    if (!f)
      errorAbort("Cannot open file: " + fileName + "\n");
  }


  inline void openOutFile(ofstream & f, string const & fileName)
  {
    f.open( fileName.c_str() );
    if (!f)
      errorAbort("Cannot open file: " + fileName + "\n");
  }



} // end namespace dgRaph

#endif  //__utils_h
