/*
 * readmatrix.cpp
 *
 *  Created on: 19-oct-2017
 *      Author: M. El-Kebir
 */

#include "readmatrix.h"

ReadMatrix::ReadMatrix()
  : BaseMatrix()
  , _var()
  , _ref()
  , _cnStates()
  , _cnStateToMu()
{
}

ReadMatrix::ReadMatrix(int m, int n)
  : BaseMatrix(m, n)
  , _var(m, IntVector(n, 0))
  , _ref(m, IntVector(n, 0))
  , _cnStates(m, CopyNumberStateMatrix(n))
  , _cnStateToMu(m, CopyNumberStateToMuMapVector(n))
{
}

ReadMatrix ReadMatrix::downSampleCharacters(int new_n) const
{
  if (new_n >= _n || new_n == -1)
  {
    return *this;
  }
  
  IntVector snvIndices(_n, 0);
  for (int i = 1; i < _n; ++i)
  {
    snvIndices[i] = snvIndices[i-1] + 1;
  }
  std::shuffle(snvIndices.begin(), snvIndices.end(), g_rng);
  snvIndices.erase(snvIndices.begin() + new_n, snvIndices.end());
  
  return downSampleCharacters(snvIndices);
}

ReadMatrix ReadMatrix::downSampleCharacters(const IntVector& snvIndices) const
{
  const int new_n = snvIndices.size();
  
  ReadMatrix newR;
  newR._m = _m;
  newR._n = new_n;
  
  newR._var = IntMatrix(_m,
                        IntVector(new_n, 0));
  newR._ref = IntMatrix(_m,
                        IntVector(new_n, 0));
  newR._cnStates = CopyNumberState3Matrix(_m,
                                          CopyNumberStateMatrix(new_n));
  newR._cnStateToMu = CopyNumberStateToMuMapMatrix(_m,
                                                   CopyNumberStateToMuMapVector(new_n));
  
  newR._indexToSample = _indexToSample;
  newR._sampleToIndex = _sampleToIndex;
  newR._indexToCharacter = StringVector(new_n);
  for (int i = 0; i < new_n; ++i)
  {
    int old_i = snvIndices[i];
    newR._indexToCharacter[i] = _indexToCharacter[old_i];
    newR._characterToIndex[newR._indexToCharacter[i]] = i;
    for (int p = 0; p < _m; ++p)
    {
      newR._var[p][i] = _var[p][old_i];
      newR._ref[p][i] = _ref[p][old_i];
      newR._cnStateToMu[p][i] = _cnStateToMu[p][old_i];
      newR._cnStates[p][i] = _cnStates[p][old_i];
    }
  }
  
  return newR;
}

int ReadMatrix::getNumberOfObservations() const
{
  // var + ref
  int res = _m * _n * 2;
  
  for (int i = 0; i < _n; ++i)
  {
    for (int p = 0; p < _m; ++p)
    {
      for (const CopyNumberState& cnState : getCopyNumberStates(p, i))
      {
        if (cnState._mu > 0)
        {
          res += 3;
        }
      }
    }
  }
  
  return res;
}

IntDoublePair ReadMatrix::computeSangerCCF(int ref,
                                           int var,
                                           double purity,
                                           const CopyNumberStateVector& cnStates)
{
  return computeSangerCCFs(ref, var, purity, cnStates)[0];
}

IntDoublePairVector ReadMatrix::computeSangerCCFs(int ref,
                                                  int var,
                                                  double purity,
                                                  const CopyNumberStateVector& cnStates)
{
  int maxCN = -1;
  double denominator = 0;
  for (const CopyNumberState& cnState : cnStates)
  {
    int cn = std::max(cnState._x, cnState._y);
    denominator += (cnState._x + cnState._y) * cnState._mu;
    if (cn > maxCN)
    {
      maxCN = cn;
    }
  }
  
  const double observedVAF = (double) var / (double) (ref + var);
  
  struct Triple {
    int _nChr;
    double _CCF;
    double _VAF;
    double _likelihood;
  };
  
  std::vector<Triple> grid;
  double maxLikelihood = -1 * std::numeric_limits<double>::max();
  for (int n_chr = 1; n_chr <= maxCN; ++n_chr)
  {
    Triple point;
    point._nChr = n_chr;
    point._VAF = observedVAF;
    point._CCF = (denominator  * observedVAF) / (n_chr * purity);
    if (point._CCF > 1)
    {
      point._CCF = 1;
      point._VAF = (n_chr * purity) / denominator;
      assert(0 <= point._VAF && point._VAF <= 1);
    }
    point._likelihood = lgamma(ref + var + 1)
                        - (lgamma(ref + 1) + lgamma(var + 1))
                        + var * log(point._VAF)
                        + ref * log(1 - point._VAF);
    maxLikelihood = std::max(point._likelihood, maxLikelihood);
    grid.push_back(point);
  }
  
  IntDoublePairVector res;
  for (const Triple& point : grid)
  {
    if (!g_tol.different(point._likelihood, maxLikelihood))
    {
      res.push_back(IntDoublePair(point._nChr, point._CCF));
    }
  }
  
  return res;
}

std::ostream& operator<<(std::ostream& out,
                         const ReadMatrix& R)
{
  out << R._m << " #samples" << std::endl;
  out << R._n << " #characters" << std::endl;
  out << "#sample_index\tsample_label\tcharacter_index\tcharacter_label\tref\tvar" << std::endl;
  for (int p = 0; p < R._m; ++p)
  {
    const std::string& pStr = R.indexToSample(p);
    for (int c = 0; c < R._n; ++c)
    {
      const std::string& cStr = R.indexToCharacter(c);
      out << p << "\t" << pStr << "\t"
          << c << "\t" << cStr << "\t"
          << R._ref[p][c] << "\t" << R._var[p][c];

      const ReadMatrix::CopyNumberStateVector& cnStates = R._cnStates[p][c];
      for (const ReadMatrix::CopyNumberState& cnState : cnStates)
      {
        out << "\t" << cnState._x
            << "\t" << cnState._y
            << "\t" << cnState._mu;
      }
      out << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in,
                         ReadMatrix& R)
{
  g_lineNumber = 0;
  std::string line;
  getline(in, line);
  
  int m = -1;
  int n = -1;
  
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: m should be nonnegative");
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: n should be nonnegative");
  }
  
  R._m = m;
  R._n = n;
  
  R._indexToSample = StringVector(m);
  R._indexToCharacter = StringVector(n);
  
  R._var = IntMatrix(m, IntVector(n, 0));
  R._ref = IntMatrix(m, IntVector(n, 0));
  R._cnStates = ReadMatrix::CopyNumberState3Matrix(m, ReadMatrix::CopyNumberStateMatrix(n));
  R._cnStateToMu = ReadMatrix::CopyNumberStateToMuMapMatrix(m, ReadMatrix::CopyNumberStateToMuMapVector(n));
  
  std::vector<std::vector<bool> > present(m, std::vector<bool>(n, false));
  while (in.good())
  {
    getline(in, line);
    if (line == "" || line[0] == '#')
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() < 6)
    {
      throw std::runtime_error(getLineNumber() + "Error: Fewer than 6 entries");
    }

    int p = -1;
    std::string pStr;
    int c = -1;
    std::string cStr;
    int ref = -1;
    int var = -1;
    try
    {
      p = boost::lexical_cast<int>(s[0]);
      pStr = s[1];
      c = boost::lexical_cast<int>(s[2]);
      cStr = s[3];
      ref = boost::lexical_cast<int>(s[4]);
      var = boost::lexical_cast<int>(s[5]);
    }
    catch (boost::bad_lexical_cast& e)
    {
      throw std::runtime_error(getLineNumber() + "Error: " + e.what());
    }
    
    if (!(0 <= p && p < m))
    {
      throw std::runtime_error(getLineNumber() + "Error: Invalid sample index '" + s[0] + "'");
    }
    
    if (!(0 <= c && c < n))
    {
      throw std::runtime_error(getLineNumber() + "Error: Invalid character index '" + s[2] + "'");
    }
    
    if (ref < 0)
    {
      throw std::runtime_error(getLineNumber() + "Error: Invalid number of reference reads '" + s[4] + "'");
    }
    
    if (var < 0)
    {
      throw std::runtime_error(getLineNumber() + "Error: Invalid number of variant reads '" + s[5] + "'");
    }
    
    if (present[p][c])
    {
      throw std::runtime_error(getLineNumber() + "Error: Duplicate character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    int offset = 6;
    int t = (s.size() - offset) / 3;
    int max_x = -1, max_y = -1;
    for (int i = 0; i < t; ++i)
    {
      ReadMatrix::CopyNumberState cnState;
      cnState._x = boost::lexical_cast<int>(s[offset + 3*i]);
      cnState._y = boost::lexical_cast<int>(s[offset + 3*i + 1]);
      cnState._mu = boost::lexical_cast<double>(s[offset + 3*i + 2]);
        
      if (cnState._x < cnState._y)
      {
        throw std::runtime_error(getLineNumber() +
                                 "Error: Found unsorted copy-number state (" + std::to_string(cnState._x) + ", " + std::to_string(cnState._y)
                                 + "); copy-number state must be ordered s.t. first copy number should always be the higher");
      }
      
      if (cnState._x > max_x) max_x = cnState._x;
      if (cnState._y > max_y) max_y = cnState._y;
      
      IntPair xy(cnState._x, cnState._y);
      
      if (R._cnStateToMu[p][c].count(xy) == 0)
      {
        R._cnStates[p][c].push_back(cnState);
        R._cnStateToMu[p][c][xy] = 0;
      }
      else
      {
        for (ReadMatrix::CopyNumberState& cnState2 : R._cnStates[p][c])
        {
          if (cnState2._x == cnState._x && cnState2._y == cnState._y)
          {
            cnState2._mu += cnState._mu;
          }
        }
      }
      
      
      R._cnStateToMu[p][c][xy] += cnState._mu;
    }
    
    present[p][c] = true;
    R._indexToSample[p] = pStr;
    R._sampleToIndex[pStr] = p;
    R._indexToCharacter[c] = cStr;
    R._characterToIndex[cStr] = c;
    R._var[p][c] = var;
    R._ref[p][c] = ref;
  
    present[p][c] = true;
  }
  
  R.consolidate();
  
  return in;
}
