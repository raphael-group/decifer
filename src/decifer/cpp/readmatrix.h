/*
 * readmatrix.h
 *
 *  Created on: 19-oct-2017
 *      Author: M. El-Kebir
 */

#ifndef READMATRIX_H
#define READMATRIX_H

#include "utils.h"
#include "basematrix.h"

class ReadMatrix : public BaseMatrix
{
public:
  struct CopyNumberState
  {
    /// Number of maternal copies
    int _x;
    /// Number of paternal copies
    int _y;
    /// Copy-number mixture proportion
    double _mu;
  };
  
  typedef std::vector<CopyNumberState> CopyNumberStateVector;
  
public:
  /// Constructor
  ReadMatrix();
  
  /// Constructor with dimensions
  ///
  /// @param m Number of samples
  /// @param n Number of characters
  ReadMatrix(int m, int n);
  
  /// Return number of observations
  int getNumberOfObservations() const;
  
  /// Downsample characters
  ///
  /// @param new_n New number of characters
  ReadMatrix downSampleCharacters(int new_n) const;
  
  /// Downsample characters
  ///
  /// @param snvIndices SNV indices to retain
  ReadMatrix downSampleCharacters(const IntVector& snvIndices) const;
  
  static IntDoublePair computeSangerCCF(int ref,
                                        int var,
                                        double purity,
                                        const CopyNumberStateVector& cnStates);
  
  static IntDoublePairVector computeSangerCCFs(int ref,
                                               int var,
                                               double purity,
                                               const CopyNumberStateVector& cnStates);
  
  /// Return whether specified character is copy neutral
  ///
  /// @param i Character
  bool isCopyNeutral(int i) const
  {
    assert(0 <= i && i < _n);
    
    for (int p = 0; p < _m; ++p)
    {
      for (const CopyNumberState& cnState : _cnStates[p][i])
      {
        if (cnState._x != 1 && cnState._y != 1 && cnState._mu > 0)
        {
          return false;
        }
      }
    }
    
    return true;
  }
  
  /// Return number of variant reads
  ///
  /// @param p Sample
  /// @param i Character
  int getVar(int p, int i) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= i && i < _n);
    
    return _var[p][i];
  }
  
  /// Return number of reference reads
  ///
  /// @param p Sample
  /// @param i Character
  int getRef(int p, int i) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= i && i < _n);
    
    return _ref[p][i];
  }
  
  /// Set number of variant reads
  ///
  /// @param p Sample
  /// @param i Character
  /// @param a_pi Number of variant reads
  void setVar(int p, int i, int a_pi)
  {
    _var[p][i] = a_pi;
  }
  
  /// Set number of reference reads
  ///
  /// @param p Sample
  /// @param i Character
  /// @param r_pi Number of reference reads
  void setRef(int p, int i, int r_pi)
  {
    _ref[p][i] = r_pi;
  }
  
  /// Set copy number mixture proportion
  ///
  /// @param p Sample
  /// @param i Character
  /// @param x Number of maternal copies
  /// @param y Number of paternal copies
  /// @param mu Copy number mixture proportion
  void setCopyNumberState(int p, int i,
                          int x, int y, double mu)
  {
    IntPair xy(x, y);
    _cnStateToMu[p][i][xy] = mu;

    CopyNumberState cnState;
    cnState._x = x;
    cnState._y = y;
    cnState._mu = mu;
    _cnStates[p][i].push_back(cnState);
  }
  
  /// Return maximum read depth
  const int getMaxReadDepth() const
  {
    int res = 0;
    for (int p = 0; p < _m; ++p)
    {
      for (int i = 0; i < _n; ++i)
      {
        res = std::max(res, _var[p][i] + _ref[p][i]);
      }
    }
    return res;
  }
  
  /// Return copy number states
  ///
  /// @param p Sample
  /// @param i Character
  const CopyNumberStateVector& getCopyNumberStates(int p, int i) const
  {
    return _cnStates[p][i];
  }
  
  /// Return maximum copy number for specified character
  ///
  /// @param i Character
  int getMaxCopies(int i) const
  {
    assert(0 <= i && i < _n);
    
    int max = 0;
    for (int p = 0; p < _m; ++p)
    {
      for (const CopyNumberState& cnState : _cnStates[p][i])
      {
        max = std::max(max, std::max(cnState._x, cnState._y));
      }
    }
    
    return max;
  }
  
  /// Identify copy number mixture proportion of
  /// provided copy number state (x,y).
  /// Return 0 if copy number state is absent.
  ///
  /// @param p Sample
  /// @param i Character
  /// @param x Number of maternal copies
  /// @param y Number of paternal copies
  double getMu(int p, int i,
               int x, int y) const
  {
    IntPair xy = IntPair(x, y);
    if (_cnStateToMu[p][i].count(xy) != 1)
    {
      return 0;
    }
    else
    {
      return _cnStateToMu[p][i].find(xy)->second;
    }
  }
  
  /// Ensure that each character has the same copy number states for each sample
  void consolidate()
  {
    for (int i = 0; i < _n; ++i)
    {
      IntPairSet X_i;
      for (int p = 0; p < _m; ++p)
      {
        for (const auto& kv : _cnStateToMu[p][i])
        {
          X_i.insert(kv.first);
        }
      }
      
      for (int p = 0; p < _m; ++p)
      {
        for (const IntPair& xy : X_i)
        {
          if (_cnStateToMu[p][i].count(xy) == 0)
          {
            _cnStateToMu[p][i][xy] = 0;
            
            CopyNumberState cnState;
            cnState._x = xy.first;
            cnState._y = xy.second;
            cnState._mu = 0;
            _cnStates[p][i].push_back(cnState);
          }
        }
      }
    }
    
    // remove (x,y) if mu = 0 in all samples
    for (int i = 0; i < _n; ++i)
    {
      IntPairSet toRemove;
      for (const CopyNumberState& cnState : _cnStates[0][i])
      {
        IntPair xy(cnState._x, cnState._y);
        bool remove = true;
        for (int p = 0; p < _m; ++p)
        {
          assert(_cnStateToMu[p][i].count(xy) == 1);
          if (_cnStateToMu[p][i][xy] > 0)
          {
            remove = false;
          }
        }
        
        if (remove)
        {
          toRemove.insert(xy);
        }
      }
      
      for (const IntPair& xy : toRemove)
      {
        for (int p = 0; p < _m; ++p)
        {
          assert(_cnStateToMu[p][i].count(xy) == 1);
          _cnStateToMu[p][i].erase(xy);
        
          for (auto it = _cnStates[p][i].begin(); it != _cnStates[p][i].end(); ++it)
          {
            if (it->_x == xy.first && it->_y == xy.second)
            {
              _cnStates[p][i].erase(it);
              break;
            }
          }
        }
      }
    }
  }
  
  IntMatrix partition() const
  {
    IntMatrix result;
    
    std::map<CopyNumberStateToMuMapVector, IntVector> map;
    for (int i = 0; i < _n; ++i)
    {
      CopyNumberStateToMuMapVector cnVector;
      for (int p = 0; p < _m; ++p)
      {
        cnVector.push_back(_cnStateToMu[p][i]);
      }
      
      map[cnVector].push_back(i);
    }
        
    for (const auto& kv : map)
    {
      result.push_back(kv.second);
    }
    
    return result;
  }
  
  double getVAF(int p, int i) const
  {
    return _var[p][i] / double(_var[p][i] + _ref[p][i]);
  }
  
protected:
  typedef std::vector<CopyNumberStateVector> CopyNumberStateMatrix;
  typedef std::vector<CopyNumberStateMatrix> CopyNumberState3Matrix;
  typedef std::map<IntPair, double> CopyNumberStateToMuMap;
  typedef std::vector<CopyNumberStateToMuMap> CopyNumberStateToMuMapVector;
  typedef std::vector<CopyNumberStateToMuMapVector> CopyNumberStateToMuMapMatrix;
  
protected:
  /// Variant reads
  IntMatrix _var;
  /// Reference reads
  IntMatrix _ref;
  /// Copy number mixture proportions
  CopyNumberState3Matrix _cnStates;
  /// Copy number mixture proportions inverse map
  CopyNumberStateToMuMapMatrix _cnStateToMu;
  
  friend std::ostream& operator<<(std::ostream& out, const ReadMatrix& R);
  friend std::istream& operator>>(std::istream& in, ReadMatrix& R);
};

std::ostream& operator<<(std::ostream& out, const ReadMatrix& R);
std::istream& operator>>(std::istream& in, ReadMatrix& R);

#endif // READMATRIX_H
