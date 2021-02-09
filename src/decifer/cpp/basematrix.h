/*
 * basematrix.h
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef BASEMATRIX_H
#define BASEMATRIX_H

#include "utils.h"

/// This class models an k * m * n matrix
class BaseMatrix
{
public:
  /// Constructor
  BaseMatrix();
  
  /// Constructor with dimensions
  ///
  /// @param m Number of samples
  /// @param n Number of characters
  BaseMatrix(int m, int n);
  
  /// Return number of samples
  int getNrSamples() const
  {
    return _m;
  }
  
  /// Return number of characters
  int getNrCharacters() const
  {
    return _n;
  }
  
  /// Return whether given label corresponds to a sample
  ///
  /// @param pStr Sample label
  bool isSample(const std::string& pStr) const
  {
    return _sampleToIndex.count(pStr) == 1;
  }
  
  /// Return whether given label corresponds to a character
  ///
  /// @param cStr Character label
  bool isCharacter(const std::string& cStr) const
  {
    return _characterToIndex.count(cStr) == 1;
  }
  
  /// Return label of given sample index
  ///
  /// @param p Sample index
  const std::string& indexToSample(int p) const
  {
    assert(0 <= p && p < _m);
    
    return _indexToSample[p];
  }
  
  /// Return index of given sample label, -1 is returned if provided
  /// sample is undefined
  ///
  /// @param pStr Sample label
  int sampleToIndex(const std::string& pStr) const
  {
    if (_sampleToIndex.count(pStr) == 1)
      return _sampleToIndex.find(pStr)->second;
    else
      return -1;
  }
  
  /// Return label of given character index
  ///
  /// @param c Character index
  const std::string& indexToCharacter(int c) const
  {
    assert(0 <= c && c < _n);
    
    return _indexToCharacter[c];
  }
  
  /// Return index of given character.
  /// If the provided does not correspond to a character, -1 is returned.
  ///
  /// @param cStr Character label
  int characterToIndex(const std::string& cStr) const
  {
    if (_characterToIndex.count(cStr) == 1)
      return _characterToIndex.find(cStr)->second;
    else
      return -1;
  }
  
protected:
  /// Number of samples
  int _m;
  /// Number of two-state characters
  int _n;
  /// Index to sample label
  StringVector _indexToSample;
  /// Sample label to index
  StringToIntMap _sampleToIndex;
  /// Index to character label
  StringVector _indexToCharacter;
  /// Character label to index
  StringToIntMap _characterToIndex;
};

#endif // BASEMATRIX_H
