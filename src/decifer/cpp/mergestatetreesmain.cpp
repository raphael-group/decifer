/*
 * mergestatetreesmain.cpp
 *
 *  Created on: 10-jun-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "stategraph.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  for (int i = 1; i < argc; ++i)
  {
    std::ifstream inS(argv[i]);
    StateGraph::readStateTrees(inS);
  }
  
  StateGraph::writeStateTrees(std::cout);
  
  return 0;
}
