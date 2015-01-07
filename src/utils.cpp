/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "utils.h"

namespace phy {

  string strip(string const & s, string const & sep)
  {
    if (s.size() == 0)
      return s;

    string::size_type begin = 0, last = s.size()-1;
    string::size_type i;
    i = begin;
    while (i < s.size()) {
      if (sep.find(s[i]) == string::npos)
	{
	  begin = i;
	  break;
	}
      i++;
    }

    i = last;
    while (i >= 0) {
      if (sep.find(s[i]) == string::npos)
	{
	  last = i;
	  break;
	}
      i--;
    }
    return s.substr(begin, last+1-begin);
  }


  vector<string> split(string const & s, string const & sep)
  {
    bool inRun = false;
    vector<string> splits;
    string::size_type begin = 0;  

    for (string::size_type i = 0; i < s.size(); i++) {
      if (not inRun) {
	if (sep.find(s[i]) == string::npos) {
	  begin = i;
	  inRun = true;
	}
      }
      else
	if (sep.find(s[i]) != string::npos) {
	  splits.push_back(s.substr(begin, i-begin));
	  inRun = false;
	}
    }
    if (inRun)
      splits.push_back(s.substr(begin, s.size() - begin));

    return splits;
  }


  istream & skipLine(istream & str, unsigned maxChar)
  {
    str.ignore(maxChar, '\n');
    return str;
  }


  istream & skipWhiteSpaceAndComments(istream & str)
  {
    char c;
    while ( str.get(c) )
      if (! isspace(c) ) {  //#include <ctype.h>
	if (c == '#') // skip rest of line
	  skipLine(str, 10000);
	else {
	  str.putback(c);
	  break;
	}
      }
    return str;
  }


  // return false if first character on line is whitespace. Skips comment lines (starting with '#')
  bool moreTags(istream & str)
  {
    // skip comment lines
    while (str.peek() == '#')
      skipLine(str);

    // check that first char is not whitespace
    if (isspace( str.peek() ) or not str.good() )
      return false;
    else
      return true;
  }


  string pathPrefix(string const & filePath)
  {
    vector<string> v = split(filePath, "/");
    if (v.size() <= 1)
      return "";
    else {
      string pathPrefix = "./";
      for (unsigned i = 0; i < v.size() - 1; i++) {
	pathPrefix += v[i] + "/";
      }
      return pathPrefix;
    }
  }

} // end namespace phy
