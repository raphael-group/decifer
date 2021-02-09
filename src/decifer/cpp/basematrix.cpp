/*
 * basematrix.cpp
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#include "basematrix.h"

BaseMatrix::BaseMatrix()
  : _m(0)
  , _n(0)
  , _indexToSample()
  , _sampleToIndex()
  , _indexToCharacter()
  , _characterToIndex()
{
}

BaseMatrix::BaseMatrix(int m, int n)
  : _m(m)
  , _n(n)
  , _indexToSample(m)
  , _sampleToIndex()
  , _indexToCharacter(n)
  , _characterToIndex()
{
  char buf[1024];
  for (int c = 0; c < n; ++c)
  {
    snprintf(buf, 1024, "SNV_%d", c + 1);
    _indexToCharacter[c] = buf;
    _characterToIndex[buf] = c;
  }
  
  for (int p = 0; p < m; ++p)
  {
    snprintf(buf, 1024, "sample_%d", p + 1);
    _indexToSample[p] = buf;
    _sampleToIndex[buf] = p;
  }
}
