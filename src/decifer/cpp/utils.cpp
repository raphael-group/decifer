/*
 * utils.cpp
 *
 *  Created on: 19-oct-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <fstream>

double g_thre(1e-6);
lemon::Tolerance<double> g_tol(g_thre);

std::mt19937 g_rng(0);

VerbosityLevel g_verbosity = VERBOSE_ESSENTIAL;

int g_lineNumber = 0;

std::string getLineNumber()
{
  char buf[1024];
  
  snprintf(buf, 1024, "Line: %d. ", g_lineNumber);
  
  return std::string(buf);
}

std::istream& getline(std::istream& is, std::string& t)
{
  ++g_lineNumber;
  
  // source: http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
  t.clear();
  
  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.
  
  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();
  
  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if(sb->sgetc() == '\n')
          sb->sbumpc();
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if(t.empty())
          is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char)c;
    }
  }
}
