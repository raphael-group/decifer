/*
 * solver.cpp
 *
 *  Created on: 12-nov-2017
 *      Author: M. El-Kebir
 */

#include "solver.h"
#include "stategraph.h"
#include <iomanip>

Solver::Solver(const ReadMatrix& R,
               int k,
               int nrSegments,
               double precisionBetaBin)
  : _R(R)
  , _k(k)
  , _nrSegments(nrSegments)
  , _precisionBetaBin(precisionBetaBin)
  , _logFactorial(1, 0)
  , _scriptT(R.getNrCharacters())
  , _x()
  , _hatG()
  , _dOverallLB()
  , _dOverallUB()
  , _dLB(R.getNrCharacters())
  , _dUB(R.getNrCharacters())
  , _denominator(R.getNrCharacters(), DoubleVector(R.getNrSamples(), 0))
  , _numerator(R.getNrCharacters())
  , _xyStar(R.getNrCharacters())
  , _logBinomCoeff(R.getNrSamples(),
                   DoubleVector(R.getNrCharacters(), 0))
  , _solD(_k, DoubleVector(R.getNrSamples(), 0))
  , _solPi(_k, 1. / _k)
  , _solZ(R.getNrCharacters(), 0)
  , _logLikelihood(-std::numeric_limits<double>::max())
{
  // init binomial coefficient
  const int maxDepth = _R.getMaxReadDepth();
  for (int i = 1; i <= maxDepth; ++i)
  {
    _logFactorial.push_back(_logFactorial.back() + log(i));
  }
}

void Solver::initPWLA()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  // compute coord
  DoubleVector coord(_nrSegments);
  coord[0] = g_tol.epsilon();
  coord[_nrSegments - 1] = 1;
  
  const double delta = 1. / (_nrSegments - 1);
  for (int alpha = 1; alpha < _nrSegments - 1; ++alpha)
  {
    coord[alpha] = delta * alpha;
  }
  
  // initialize segments
  const double infeasible_log = -1e300;
  _hatG = Double5Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int scriptT_i_size = _scriptT[i].size();
    _hatG[i] = Double4Matrix(scriptT_i_size);
    for (int t = 0; t < scriptT_i_size; ++t)
    {
      _hatG[i][t] = DoubleTensor(_k);
      for (int j = 0; j < _k; ++j)
      {
        _hatG[i][t][j] = DoubleMatrix(m);
        for (int p = 0; p < m; ++p)
        {
          _hatG[i][t][j][p] = DoubleVector(_nrSegments, 0);
          
          const int var_ip = _R.getVar(p, i);
          const int ref_ip = _R.getRef(p, i);
          
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            double h = (coord[alpha] - _numerator[i][t][p]) / _denominator[i][p];
            if (!g_tol.nonZero(h))
            {
              h = 0;
            }
            else if (!g_tol.different(h, 1))
            {
              h = 1;
            }
//            if (var_ip > 0 && (h <= 0 || !g_tol.nonZero(h)))
//            {
//              _hatG[i][t][j][p][alpha] = infeasible_log;//log(g_tol.epsilon());
//            }
//            else if (ref_ip > 0 && (h >= 1 || g_tol.less(1, h)))
//            {
//              _hatG[i][t][j][p][alpha] = infeasible_log;//log(g_tol.epsilon());
//            }
//            else
            {
              double likelihood = getLogLikelihood(var_ip, ref_ip, h);
//              assert(!isinf(likelihood));
//              assert(!isnan(likelihood));
              if (likelihood < infeasible_log || isnan(likelihood) || isinf(likelihood))
              {
                _hatG[i][t][j][p][alpha] = infeasible_log;//log(g_tol.epsilon());
              }
              else
              {
                _hatG[i][t][j][p][alpha] = likelihood;
              }
            }
          }
        }
      }
    }
  }
}

StateTree Solver::convertToStateTreeFromSNVF(const ReadMatrix& R,
                                             const StateEdgeSet& T_it,
                                             const DoubleVector& f_i,
                                             const int i)
{
  assert(R.getNrSamples() == f_i.size());
  
  IntVector pi;
  std::map<StateGraph::CnaTriple, int> vertices;
  
  StateGraph::CnaTripleSet verticesPreMut, verticesPostMut;
  StateGraph::CnaTriple mutationVertex;
  StateGraph::partition(T_it, verticesPreMut, verticesPostMut, mutationVertex);
  
  assert(mutationVertex._x != -1);
  assert(mutationVertex._y != -1);
  
  vertices[StateGraph::CnaTriple(1, 1, 0)] = 0;
  pi.push_back(-1);
  for (const StateGraph::StateEdge& edge : T_it)
  {
    if (vertices.count(edge.first) == 0)
    {
      vertices[edge.first] = pi.size();
      pi.push_back(-1);
    }
    if (vertices.count(edge.second) == 0)
    {
      vertices[edge.second] = pi.size();
      pi.push_back(vertices[edge.first]);
    }
    pi[vertices[edge.second]] = vertices[edge.first];
  }
  
  StateTree S_it(pi);
  
  for (int p = 0; p < f_i.size(); ++p)
  {
    double muTotal = 0;
    for (const ReadMatrix::CopyNumberState& cnState : R.getCopyNumberStates(p, i))
    {
      muTotal += (cnState._x + cnState._y) * cnState._mu;
    }
    
    double muMut = 0;
    for (const auto& kv : vertices)
    {
      if (kv.first._x != mutationVertex._x || kv.first._y != mutationVertex._y)
      {
        muMut += kv.first._z * R.getMu(p, i, kv.first._x, kv.first._y);
      }
    }
    
    double muStar = R.getMu(p, i, mutationVertex._x, mutationVertex._y);
    
    const double h_ip = f_i[p];
    
    for (const auto& kv : vertices)
    {
      if (p == 0)
      {
        S_it.setCnState(kv.second, kv.first._x, kv.first._y, kv.first._z);
      }
      
      if (kv.first._x != mutationVertex._x || kv.first._y != mutationVertex._y)
      {
        S_it.setMixtureProportion(kv.second, R.getMu(p, i, kv.first._x, kv.first._y));
      }
      else if (kv.first._x == mutationVertex._x && kv.first._y == mutationVertex._y && kv.first._z == 1)
      {
        double s_mut = h_ip * muTotal - muMut;
        
        S_it.setMixtureProportion(kv.second, s_mut);
      }
      else
      {
        assert(kv.first._x == mutationVertex._x && kv.first._y == mutationVertex._y && kv.first._z == 0);
        
        double s_nonMut = muStar - (h_ip * muTotal - muMut);
        
        S_it.setMixtureProportion(kv.second, s_nonMut);
      }
    }
  }
  
  return S_it;
}

StateTree Solver::convertToStateTreeFromDCF(const ReadMatrix& R,
                                            const StateEdgeSet& T_it,
                                            const DoubleVector& d_i,
                                            const int i)
{
  assert(R.getNrSamples() == d_i.size());

  // parent vector
  IntVector pi;
  std::map<StateGraph::CnaTriple, int> vertices;
  
  StateGraph::CnaTripleSet verticesPreMut, verticesPostMut;
  StateGraph::CnaTriple mutationVertex;
  StateGraph::partition(T_it, verticesPreMut, verticesPostMut, mutationVertex);
  
  assert(mutationVertex._x != -1);
  assert(mutationVertex._y != -1);
  
  vertices[StateGraph::CnaTriple(1, 1, 0)] = 0;
  pi.push_back(-1);
  for (const StateGraph::StateEdge& edge : T_it)
  {
    if (vertices.count(edge.first) == 0)
    {
      vertices[edge.first] = pi.size();
      pi.push_back(-1);
    }
    if (vertices.count(edge.second) == 0)
    {
      vertices[edge.second] = pi.size();
      pi.push_back(vertices[edge.first]);
    }
    pi[vertices[edge.second]] = vertices[edge.first];
  }
  
  StateTree S_it(pi);
  for (int p = 0; p < d_i.size(); ++p)
  {
    for (const auto& kv : vertices)
    {
      if (p == 0)
      {
        S_it.setCnState(kv.second, kv.first._x, kv.first._y, kv.first._z);
      }
    }
    
    double massPreMut = 0, massPostMut = 0, massMut = 0;
    
    for (int idx = 0; idx < S_it.numVertices(); ++idx)
    {
      int x = S_it.x(idx);
      int y = S_it.y(idx);
      int z = S_it.z(idx);
      if (S_it.isPresent(idx))
      {
        if (x == mutationVertex._x && y == mutationVertex._y && z == mutationVertex._z)
        {
          massMut = R.getMu(p, i, x, y);
        }
        else if (verticesPreMut.count(StateGraph::CnaTriple(x, y, z)) == 1)
        {
          massPreMut += R.getMu(p, i, x, y);
        }
        else
        {
          assert(verticesPostMut.count(StateGraph::CnaTriple(x, y, z)) == 1);
          massPostMut += R.getMu(p, i, x, y);
        }
      }
    }
    
    double muTotal = 0;
    for (const ReadMatrix::CopyNumberState& cnState : R.getCopyNumberStates(p, i))
    {
      muTotal += (cnState._x + cnState._y) * cnState._mu;
    }
    
    // add all mu's that are not (x*,y*)
    double muMut = 0;
    for (const auto& kv : vertices)
    {
      if (kv.first._x != mutationVertex._x || kv.first._y != mutationVertex._y)
      {
        muMut += kv.first._z * R.getMu(p, i, kv.first._x, kv.first._y);
      }
    }
    
    const double d_ip = d_i[p];    
    for (const auto& kv : vertices)
    {
      if (p == 0)
      {
        S_it.setCnState(kv.second, kv.first._x, kv.first._y, kv.first._z);
      }
      
      if (kv.first._x != mutationVertex._x || kv.first._y != mutationVertex._y)
      {
        S_it.setMixtureProportion(kv.second, R.getMu(p, i, kv.first._x, kv.first._y));
      }
      else if (kv.first._x == mutationVertex._x && kv.first._y == mutationVertex._y && kv.first._z == 1)
      {
        double s_mut = d_ip - massPostMut;
        
        S_it.setMixtureProportion(kv.second, s_mut);
      }
      else
      {
        assert(kv.first._x == mutationVertex._x && kv.first._y == mutationVertex._y && kv.first._z == 0);
        
        double s_nonMut = massMut - (d_ip - massPostMut);
        
        S_it.setMixtureProportion(kv.second, s_nonMut);
      }
    }
  }
  
  return S_it;
}

void Solver::init(int i)
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  assert(0 <= i && i < n);
  
  // 0. obtain L
  const ReadMatrix::CopyNumberStateVector& cnStates = _R.getCopyNumberStates(0, i);
  // + 1 is done on purpose, allowing the following state tree for instance
  // (1,1,0) -> (2,2,0) -> (2,2,1) -> (2,0,1)
  const int maxCopies = _R.getMaxCopies(i) + 1;
  
  IntPairSet L;
  for (const ReadMatrix::CopyNumberState& cnState : cnStates)
  {
    L.insert(IntPair(cnState._x, cnState._y));
  }
  
  // 1. enumerate all state trees
  const auto& res = StateGraph::getStateTrees(L, maxCopies);
  _scriptT[i] = StateEdgeSetVector(res.begin(), res.end());
  
  // 2. initialize D
  for (int p = 0; p < m; ++p)
  {
    const ReadMatrix::CopyNumberStateVector& cnStates_p = _R.getCopyNumberStates(p, i);
    for (const ReadMatrix::CopyNumberState& cnState : cnStates_p)
    {
      _denominator[i][p] += (cnState._x + cnState._y) * cnState._mu;
    }
  }
  
  // 3. initialize C, _xyStar, _gamma, f_LB and _fUB
  _numerator[i] = DoubleMatrix(_scriptT[i].size(),
                       DoubleVector(m, 0));
  _xyStar[i] = IntPairVector(_scriptT[i].size(),
                             IntPair(-1, -1));
  _dLB[i] = DoubleMatrix(_scriptT[i].size(),
                         DoubleVector(m, 0));
  _dUB[i] = DoubleMatrix(_scriptT[i].size(),
                         DoubleVector(m, 1));
  
  for (int t = 0; t < _scriptT[i].size(); ++t)
  {
    StateGraph::StateEdgeSet& T_it = _scriptT[i][t];
    StateGraph::CnaTripleSet verticesPreMut, verticesPostMut;
    StateGraph::CnaTriple mutationVertex;
    StateGraph::partition(T_it, verticesPreMut, verticesPostMut, mutationVertex);
    _xyStar[i][t].first = mutationVertex._x;
    _xyStar[i][t].second = mutationVertex._y;
    
    for (int p = 0; p < m; ++p)
    {
      for (const auto& xyz : verticesPostMut)
      {
        if (xyz == mutationVertex) continue;
        
//        switch (_statType)
//        {
//          case CLUSTER_DCF:
            _dLB[i][t][p] += _R.getMu(p, i, xyz._x, xyz._y);
            _numerator[i][t][p] += (1 - xyz._z) * _R.getMu(p, i, xyz._x, xyz._y);
//            break;
//          case CLUSTER_CCF:
//            if (xyz._z > 0)
//            {
              _dLB[i][t][p] += _R.getMu(p, i, xyz._x, xyz._y);
              _numerator[i][t][p] += (1 - xyz._z) * _R.getMu(p, i, xyz._x, xyz._y);
//            }
//            break;
//        }
        
      }
      _dUB[i][t][p] = _dLB[i][t][p] + _R.getMu(p, i, mutationVertex._x, mutationVertex._y);
      if (!g_tol.different(_dLB[i][t][p], _dUB[i][t][p]))
      {
        _dUB[i][t][p] = _dLB[i][t][p];
      }
      
//      const double h_lb = (_dLB[i][t][p] - _numerator[i][t][p]) / _denominator[i][p];
//      if ((h_lb <= 0 || !g_tol.nonZero(h_lb)) && _R.getVar(p, i) != 0)
//      {
//        _dLB[i][t][p] += 1e-2;
//      }
//      const double h_ub = (_dUB[i][t][p] - _numerator[i][t][p]) / _denominator[i][p];
//      if ((h_ub >= 1 || g_tol.less(1, h_ub)) && _R.getRef(p, i) != 0)
//      {
//        _dUB[i][t][p] -= 1e-2;
//      }
    }
  }
  
  // 4. initialize _logBinomCoeff
  for (int p = 0; p < m; ++p)
  {
    int var_pi = _R.getVar(p, i);
    int ref_pi = _R.getRef(p, i);
    int tot_pi = var_pi + ref_pi;
    _logBinomCoeff[p][i] = _logFactorial[tot_pi] - (_logFactorial[var_pi] + _logFactorial[ref_pi]);
  }
}

void Solver::init()
{
  const int n = _R.getNrCharacters();
  
  for (int i = 0; i < n; ++i)
  {
    init(i);
  }
  
  return;
}

double Solver::getLogLikelihood(int var, int ref, double f) const
{
  if (!g_tol.nonZero(f))
  {
    f = 1e-3;
  }
  if (!g_tol.different(1, f))
  {
    f = 1 - 1e-3;
  }

  if (!g_tol.nonZero(f) && var == 0)
  {
    return 0;
  }
  if (!g_tol.different(1, f) && ref == 0)
  {
    return 0;
  }
  
  if (_precisionBetaBin > 0)
  {
    double val = lgamma(var + ref + 1) + lgamma(var + _precisionBetaBin * f) + lgamma(ref + _precisionBetaBin * (1-f)) + lgamma(_precisionBetaBin);
    val -= (lgamma(var + 1) + lgamma(ref + 1) + lgamma(var + ref + _precisionBetaBin) + lgamma(_precisionBetaBin*f) + lgamma(_precisionBetaBin*(1-f)));
    
    return val;
  }
  else
  {
    return lgamma(var + ref + 1) - lgamma(var + 1) - lgamma(ref + 1) + var * log(f) + ref * log(1 - f);
  }
}

double Solver::getLogLikelihood(int i) const
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  assert(0 <= i && i < n);
  
  const int maxCopyNumber = _R.getMaxCopies(i);
  
  double sum_i = 0;
  const StateEdgeSetVector& scriptT_i = _scriptT[i];
  const int size_scriptT_i = scriptT_i.size();
  for (int t = 0; t < size_scriptT_i; ++t)
  {
    const StateEdgeSet& T_it = scriptT_i[t];
    for (int j = 0; j < _k; ++j)
    {
      if (_solPi[j] == 0)
        continue;
      
      if (!isEnabled(i, t, j))
        continue;
      
      double prod = _solPi[j] / size_scriptT_i;
      double log_prod = log(prod);
      for (int p = 0; p < m; ++p)
      {
        const ReadMatrix::CopyNumberStateVector& cnStates_pi = _R.getCopyNumberStates(p, i);
        const int var_pi = _R.getVar(p, i);
        const int ref_pi = _R.getRef(p, i);
        
        if (isFeasible(_solD[j][p],
                       _xyStar[i][t],
                       cnStates_pi,
                       T_it,
                       maxCopyNumber))
        {
          const double h = (_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p];
//          assert(_denominator[i][p] != 0);
//          if (!((!(h <= 0 || !g_tol.nonZero(h)) || var_pi == 0) && (!(h >= 1 || g_tol.less(1, h)) || ref_pi == 0)))
//          {
//            prod = 0;
//            break;
//          }
//          else
          {
            const double val = getLogLikelihood(var_pi, ref_pi, h);
            log_prod += val;
//            prod *= exp(val);
//            prod = std::max(std::numeric_limits<double>::min(), prod);
          }
        }
        else
        {
          prod = 0;
          break;
        }
      }
      if (prod != 0)
      {
        //sum_i += prod;
        prod = exp(log_prod);
        prod = std::max(std::numeric_limits<double>::min(), prod);
        sum_i += prod;
      }
    }
  }
  return log(sum_i);
}

double Solver::updateLogLikelihood()
{
  const int n = _R.getNrCharacters();
  
  double oldLogLikelihood = _logLikelihood;
  
  _logLikelihood = 0;
  for (int i = 0; i < n; ++i)
  {
    double L_i = getLogLikelihood(i);
    _logLikelihood += L_i;
//    std::cout << i << ": " << L_i << std::endl;
  }
//  std::cout << std::endl;
  
  return _logLikelihood - oldLogLikelihood;
}

bool Solver::isFeasible(const double f_pj,
                        const IntPair& xyStar,
                        const ReadMatrix::CopyNumberStateVector& cnStates,
                        const StateEdgeSet& T_it,
                        const int maxCopyNumber) const
{
  DoubleTensor s(maxCopyNumber + 1,
                 DoubleMatrix(maxCopyNumber + 1,
                              DoubleVector(maxCopyNumber + 1, 0)));
  
  // 2a. get vertices of state tree
  StateGraph::CnaTripleSet vertices;
  vertices.insert(StateGraph::CnaTriple(1,1,0));
  
  StateGraph::StateEdge mutEdge;
  for (const StateGraph::StateEdge& st : T_it)
  {
    vertices.insert(st.first);
    vertices.insert(st.second);
    
    // check whether st is a mutation edge
    if (st.second._x == xyStar.first && st.second._y == xyStar.second && st.second._x == st.first._x && st.second._y == st.first._y)
    {
      assert(st.second._z > 0);
      mutEdge = st;
    }
  }
  
  StateGraph::CnaTripleSet verticesPreMutS, verticesPostMutS;
  StateGraph::CnaTriple mutationVertex;
  StateGraph::partition(T_it, verticesPreMutS, verticesPostMutS, mutationVertex);
  
  // 2b. compute mass pre mutation, and mass post mut
  double massPreMut = 0, massPostMut = 0, massMut = 0;
  for (const ReadMatrix::CopyNumberState& cnState : cnStates)
  {
    if (cnState._x == xyStar.first && cnState._y == xyStar.second)
    {
      massMut = cnState._mu;
    }
  }
  
  for (const StateGraph::CnaTriple& xyz : vertices)
  {
    if (xyz != mutationVertex)
    {
      for (const ReadMatrix::CopyNumberState& cnState : cnStates)
      {
        if (cnState._x == xyz._x && cnState._y == xyz._y)
        {
          if (verticesPostMutS.count(xyz) == 1)
          {
//            if (_statType == CLUSTER_DCF || xyz._z > 0)
            {
              massPostMut += cnState._mu;
            }
            s[xyz._x][xyz._y][xyz._z] = massPostMut;
          }
          else
          {
            assert(verticesPreMutS.count(xyz) == 1);
            assert(xyz._z == 0);
            massPreMut += cnState._mu;
            s[xyz._x][xyz._y][xyz._z] = massPreMut;
          }
        }
      }
    }
  }
  
  if (g_tol.less(f_pj, massPostMut))
  {
    return false;
  }
  
  double slack = f_pj - massPostMut;
  if (g_tol.less(massMut, slack))
  {
    return false;
  }
  
  s[mutEdge.second._x][mutEdge.second._y][mutEdge.second._z] = slack;
  s[mutEdge.first._x][mutEdge.first._y][mutEdge.first._z] = massMut - slack;
  
  return true;
}

void Solver::writeSolution(std::ostream& out) const
{
  const int m = _R.getNrSamples();
  
  out << _k << " #clusters" << std::endl;
  out << m << " #samples" << std::endl;
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      if (p != 0)
      {
        out << " ";
      }
      out << _solD[j][p];
    }
    out << std::endl;
  }
  
  for (int j = 0; j < _k; ++j)
  {
    out << _solPi[j] << std::endl;
  }
}

void Solver::readSolution(std::istream& in,
                          DoubleMatrix& outF, DoubleVector& outPi)
{
  g_lineNumber = 0;
  int k = -1;
  int m = -1;
  
  g_lineNumber = 0;
  std::string line;
  getline(in, line);
  
  std::stringstream ss(line);
  ss >> k;
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> m;
  
  DoubleMatrix f(k, DoubleVector(m, -1));
  for (int j = 0; j < k; ++j)
  {
    getline(in, line);
    ss.clear();
    ss.str(line);
    for (int p = 0; p < m; ++p)
    {
      ss >> f[j][p];
      if (!(0 <= f[j][p] && f[j][p] <= 1))
      {
        throw std::runtime_error(getLineNumber() + "Error: Invalid DCF");
      }
    }
  }
  
  DoubleVector pi(k, -1);
  for (int j = 0; j < k; ++j)
  {
    getline(in, line);
    ss.clear();
    ss.str(line);
    
    ss >> pi[j];
    if (!(0 <= pi[j] && pi[j] <= 1))
    {
      throw std::runtime_error(getLineNumber() + "Error: Invalid pi.");
    }
  }
  
  outF = f;
  outPi = pi;
}

bool Solver::readSolution(std::istream& in)
{
  DoubleMatrix f;
  DoubleVector pi;
  try
  {
    readSolution(in, f, pi);
    if (f.size() != _k)
    {
      std::cerr << "Error: Invalid k" << std::endl;
      return false;
    }
    if (f[0].size() != _R.getNrSamples())
    {
      std::cerr << "Error: Invalid m" << std::endl;
      return false;
    }
    _solD = f;
    _solPi = pi;
    updateLogLikelihood();
    return true;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
}

void Solver::writeClustering(std::ostream& out) const
{
  const int n = _R.getNrCharacters();
  
  for (int i = 0; i < n; ++i)
  {
    out << i << "\t" << _R.indexToCharacter(i) << "\t" << _solZ[i] << std::endl;
  }
}

void Solver::writeMutationProperties(std::ostream& out) const
{
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  
  out << "SNV";
  out << "\tcluster";
  
  for (int p = 0; p < m; ++p)
  {
//    if (_statType == CLUSTER_DCF)
    {
      out << "\tdcf" << p;
    }
//    else
//    {
//      out << "\tccf" << p;
//    }
  }
  out << std::endl;
  
  for (int i = 0; i < n; ++i)
  {
    out << i;
    out << "\t" << _solZ[i];
    
    for (int p = 0; p < m; ++p)
    {
      out << "\t" << std::setprecision(2) << _solD[_solZ[i]][p];
    }
    out << std::endl;
  }
}
