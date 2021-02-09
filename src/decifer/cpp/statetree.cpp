/*
 *  statetree.cpp
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#include "statetree.h"

#include <lemon/connectivity.h>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iomanip>

StateTree::StateTree()
  : _k(0)
  , _S()
  , _root(lemon::INVALID)
  , _label(_S)
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
  , _cnState(_S)
{
}

StateTree::StateTree(int k)
  : _k(k)
  , _S()
  , _root(lemon::INVALID)
  , _label(_S)
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
  , _cnState(_S)
{
  IntVector pi(_k, -1);
  for (int i = 1; i < _k; ++i)
  {
    pi[i] = i - 1;
  }
  
  init(pi);
}

void StateTree::writeSummary(const ReadMatrix& R,
                             const StateTreeVector& sol,
                             std::ostream& out)
{
  const int n = R.getNrCharacters();
  const int m = R.getNrSamples();
  
  out << "SNV_index\t"
      << "SNV_label\t"
      << "sample_index\t"
      << "sample_label\t"
      << "ref\t"
      << "alt\t"
      << "VAF\t"
      << "sanger_n_mut\t"
      << "sanger_CCF\t"
      << "decifer_CCF\t"
      << "decifer_DCF\t"
      << "decifer_cluster_CCF\t"
      << "decifer_cluster_DCF\t"
//      << "decifer_cluster_index\t"
//      << "decifer_state_tree_index\t"
      << "nr_nonzero_z"
      << std::endl;

  for (int i = 0; i < n; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      int ref = R.getRef(p, i);
      int var = R.getVar(p, i);
      double vaf = R.getVAF(p, i);
      double purity = 1;
      auto sanger = R.computeSangerCCF(ref, var,
                                       purity, R.getCopyNumberStates(p, i));
      int sanger_n_mut = sanger.first;
      double sangerCCF = sanger.second;
      
      IntSet zSet;
      IntPairSet mutTumorSet;
      IntPairSet nonMutTumorSet;
      const StateTree& T_i = sol[i];
      for (int j = 0; j != T_i.numVertices(); ++j)
      {
        if (T_i.isPresent(j))
        {
          int x_j = T_i.x(j);
          int y_j = T_i.y(j);
          int z_j = T_i.z(j);
          if (z_j > 0)
          {
            zSet.insert(z_j);
          }
          if (x_j != 1 || y_j != 1)
          {
            if (z_j > 0)
            {
              mutTumorSet.insert(IntPair(x_j, y_j));
            }
            else
            {
              nonMutTumorSet.insert(IntPair(x_j, y_j));
            }
          }
        }
      }
      
      double deciferCCF = T_i.maxLikelihoodCCF(p, vaf);
      double deciferDCF = T_i.maxLikelihoodDCF(p, vaf);
      double deciferClusterDCF = T_i.dcf(p);
      double deciferClusterCCF = T_i.cf(p);
      
      out << i << "\t" << R.indexToCharacter(i) << "\t"
          << p << "\t" << R.indexToSample(p) << "\t"
          << ref << "\t" << var << "\t"
          << R.getVAF(p, i) << "\t"
          << sanger_n_mut << "\t" << sangerCCF << "\t"
          << deciferCCF << "\t" << deciferDCF << "\t"
          << deciferClusterCCF << "\t" << deciferClusterDCF << "\t"
          << zSet.size()
          << std::endl;
    }
  }
}

Node StateTree::computeMaxLikelihoodStateProportions(DoubleNodeMap& s,
                                                     int p,
                                                     double vaf) const
{
  // 1. identify mutation vertex
  int star = -1;
  int x_star = -1;
  int y_star = -1;
  double mu_star = 0;
  Node v_star = lemon::INVALID;
  
  for (NodeIt v_j(_S); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (isPresent(j))
    {
      int i = parent(j);
      if (i == -1 || !isPresent(i)) continue;
      if (z(i) == 0 && z(j) == 1)
      {
        star = j;
        x_star = x(i);
        y_star = y(i);
        v_star = v_j;
        mu_star = _cnState[_stateToNode[i]]._s[p] + _cnState[v_j]._s[p];
        break;
      }
    }
  }
  
  assert(star != -1);
  
  double C = mu_star * (x_star + y_star);
  double D = 0;
  for (NodeIt v_j(_S); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (isPresent(j) && (x(j) != x_star || y(j) != y_star))
    {
      C += _cnState[v_j]._s[p] * (x(j) + y(j));
      D += _cnState[v_j]._s[p] * z(j);
    }
  }
  
  for (NodeIt v_j(_S); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (x_star == x(j) && y_star == y(j))
    {
      if (j == star)
      {
        s[v_j] = std::min(mu_star, vaf * C - D);
      }
      else
      {
        s[v_j] = std::max(0., mu_star - vaf * C + D);
      }
    }
    else
    {
      s[v_j] = _cnState[v_j]._s[p];
    }
  }
  
  return v_star;
}

double StateTree::maxLikelihoodCCF(int p,
                                   double vaf) const
{
  DoubleNodeMap s(_S, 0.);
  computeMaxLikelihoodStateProportions(s, p, vaf);
  
  double ccf = 0;
  for (NodeIt v_j(_S); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (isPresent(j) && z(j) > 0)
    {
      ccf += s[v_j];
    }
  }
  
  return ccf;
}

double StateTree::maxLikelihoodDCF(int p,
                                   double vaf) const
{
  DoubleNodeMap s(_S, 0.);
  Node v_star = computeMaxLikelihoodStateProportions(s, p, vaf);
  int star = _nodeToState[v_star];
  
  double dcf = 0;
  for (NodeIt v_j(_S); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
//    if (isPresent(j))
//      std::cout << "s_{" << x(j) << "," << y(j) << "," << z(j) << "} = " << s[v_j] << " " << _cnState[v_j]._s[p] << std::endl;
    if (isPresent(j) && isDescendant(j, star))
    {
      dcf += s[v_j];
    }
  }
  
  return dcf;
}

void StateTree::init(const IntVector& pi)
{
  // state 0 should be the root
  assert(pi[0] == -1);
  char buf[1024];
  
  _root = _S.addNode();
  _stateToNode[0] = _root;
  _nodeToState[_root] = 0;
  snprintf(buf, 1024, "%d", 0);
  _label[_root] = buf;
  
  // init nodes of S_c
  for (int i = 1; i < _k; ++i)
  {
//    if (0 <= pi[i] && pi[i] < _k)
//    if (pi[i] != -2)
    {
      Node v_i = _S.addNode();
      _stateToNode[i] = v_i;
      _nodeToState[v_i] = i;
      snprintf(buf, 1024, "%d", i);
      _label[v_i] = buf;
    }
  }
  
  // init edges of S_c
  for (int i = 1; i < _k; ++i)
  {
    int pi_i = pi[i];
    
    //assert(0 <= pi_i < _k);
    if (0 <= pi_i && pi_i < _k)
    {
      _S.addArc(_stateToNode[pi_i], _stateToNode[i]);
    }
  }
  
  initD(_root);
  
  assert(lemon::dag(_S));
}
  
StateTree::StateTree(const IntVector& pi)
  : _k(pi.size())
  , _S()
  , _root(lemon::INVALID)
  , _label(_S, "")
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
  , _cnState(_S)
{
  init(pi);
}
  
StateTree::StateTree(const StateTree& other)
  : _k(other._k)
  , _S()
  , _root(lemon::INVALID)
  , _label(_S)
  , _stateToNode(_k, lemon::INVALID)
  , _nodeToState(_S)
  , _D(_S)
  , _cnState(_S)
{
  lemon::digraphCopy(other._S, _S)
    .node(other._root, _root)
    .nodeMap(other._nodeToState, _nodeToState)
    .nodeMap(other._D, _D)
    .nodeMap(other._label, _label)
    .nodeMap(other._cnState, _cnState)
    .run();
  
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    _stateToNode[_nodeToState[v_i]] = v_i;
  }
}
  
StateTree& StateTree::operator=(const StateTree& other)
{
  if (this != &other)
  {
    _k = other._k;
    _stateToNode = NodeVector(_k, lemon::INVALID);
    
    lemon::digraphCopy(other._S, _S)
      .node(other._root, _root)
      .nodeMap(other._nodeToState, _nodeToState)
      .nodeMap(other._label, _label)
      .nodeMap(other._D, _D)
      .nodeMap(other._cnState, _cnState)
      .run();
    
    for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
    {
      _stateToNode[_nodeToState[v_i]] = v_i;
    }
  }
  return *this;
}
  
void StateTree::writeEdgeList(std::ostream& out) const
{
  bool first = true;
  for (ArcIt a_ij(_S); a_ij != lemon::INVALID; ++a_ij)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      out << " ; ";
    }
    
    Node v_i = _S.source(a_ij);
    Node v_j = _S.target(a_ij);
    
    out << _label[v_i] << " -> " << _label[v_j];
  }
}

void StateTree::writeDOT(const int m,
                         std::ostream& out) const
{
  out << "digraph S {" << std::endl;
  out << "\tlabel=\"VAF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << vaf(p);
  }
  out << "\\nCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << cf(p);
  }
  out << "\\nDCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << dcf(p);
  }
  out << "\"" << std::endl;
  
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;
    
    out << "\t" << i << " [label=\"("
        << _cnState[v_i]._x << ","
        << _cnState[v_i]._y << ","
        << _cnState[v_i]._z << ")";
    for (double s : _cnState[v_i]._s)
    {
      out << std::setprecision(2) << "\\n" << s;
    }
    out << "\"]" << std::endl;
  }
  
  for (ArcIt a_ij(_S); a_ij != lemon::INVALID; ++a_ij)
  {
    int i = _nodeToState[_S.source(a_ij)];
    int j = _nodeToState[_S.target(a_ij)];
    
    out << "\t" << i << " -> " << j << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void StateTree::writeDOT(std::ostream& out) const
{
  out << "digraph S {" << std::endl;
  
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;
    
    out << "\t" << i << " [label=\"" << _label[v_i];
    for (double s : _cnState[v_i]._s)
    {
      out << std::setprecision(2) << "\\n" << s;
    }
    out << "\"]" << std::endl;
  }
  
  for (ArcIt a_ij(_S); a_ij != lemon::INVALID; ++a_ij)
  {
    int i = _nodeToState[_S.source(a_ij)];
    int j = _nodeToState[_S.target(a_ij)];
    
    out << "\t" << i << " -> " << j << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void StateTree::initD(Node v_i)
{
  IntSet& D_i = _D[v_i];
  D_i.clear();
  
  D_i.insert(_nodeToState[v_i]);
  
  for (OutArcIt a(_S, v_i); a != lemon::INVALID; ++a)
  {
    Node v_j = _S.target(a);
    
    initD(v_j);
    D_i.insert(_D[v_j].begin(), _D[v_j].end());
  }
}

void StateTree::sampleMixtureProportions(double f, double purity)
{
  const int mutState = mutationState();
  const IntSet& D_mutState = D(mutState);
  std::gamma_distribution<> gamma_dist(1, 1);
  
  NodeSet postMutNodes;
  NodeSet preMutNodes;
  double sum_pre = 0;
  double sum_post = 0;
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    const int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;

    _cnState[v_i]._s.push_back(gamma_dist(g_rng));
    if (D_mutState.count(i) == 1)
    {
      postMutNodes.insert(v_i);
      sum_post += _cnState[v_i]._s.back();
    }
    else
    {
      preMutNodes.insert(v_i);
      sum_pre += _cnState[v_i]._s.back();
    }
  }
  
  // zero out nodes with small weights
  for (Node v_i : postMutNodes)
  {
    if ((_cnState[v_i]._s.back() / sum_post) < 0.05)
    {
      sum_post -= _cnState[v_i]._s.back();
      _cnState[v_i]._s.back() = 0;
    }
  }
  for (Node v_i : preMutNodes)
  {
    if ((_cnState[v_i]._s.back() / sum_pre) < 0.05)
    {
      sum_pre -= _cnState[v_i]._s.back();
      _cnState[v_i]._s.back() = 0;
    }
  }
  
  for (Node v_i : postMutNodes)
  {
    _cnState[v_i]._s.back() /= sum_post;
    _cnState[v_i]._s.back() *= f;
    if (!g_tol.nonZero(_cnState[v_i]._s.back()))
    {
      _cnState[v_i]._s.back() = 0;
    }
  }
  
  for (Node v_i : preMutNodes)
  {
    _cnState[v_i]._s.back() /= sum_pre;
    _cnState[v_i]._s.back() *= (1 - f - (1 - purity));
    if (!g_tol.nonZero(_cnState[v_i]._s.back()))
    {
      _cnState[v_i]._s.back() = 0;
    }
  }
  
  _cnState[_root]._s.back() += 1 - purity;
  if (!g_tol.nonZero(_cnState[_root]._s.back()))
  {
    _cnState[_root]._s.back() = 0;
  }
}

void StateTree::writeClusterDOT(double gamma,
                                int t,
                                std::ostream& out) const
{
  const int m = getNrSamples();
  
  out << "\tsubgraph cluster_" << t << " {" << std::endl;
  out << "\t\tlabel=\"gamma = " << gamma;
  
  out << "\\nVAF = ";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << vaf(p);
  }
  out << "\\nCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << cf(p);
  }
  out << "\\nDCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << dcf(p);
  }
  out << "\"" << std::endl;
  
  for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;
    
    out << "\t\t" << i + t * 100 << " [label=\"("
        << _cnState[v_i]._x << ","
        << _cnState[v_i]._y << ","
        << _cnState[v_i]._z << ")";
    for (double s : _cnState[v_i]._s)
    {
      out << std::setprecision(2) << "\\n" << s;
    }
    out << "\"]" << std::endl;
  }
  
  for (ArcIt a_ij(_S); a_ij != lemon::INVALID; ++a_ij)
  {
    int i = _nodeToState[_S.source(a_ij)];
    int j = _nodeToState[_S.target(a_ij)];
    
    out << "\t\t" << i + t * 100 << " -> " << j + t * 100 << std::endl;
  }
  
  out << "\t}" << std::endl;
}

std::ostream& operator<<(std::ostream& out, const StateTree& S)
{
  // output pi
  out << -1;
  for (int i = 1; i < S.k(); ++i)
  {
    out << " " << S.parent(i);
  }
  out << std::endl;
  
  // output vertices
  for (int i = 0; i < S.k(); ++i)
  {
    Node v_i = S.node(i);
    out << S._cnState[v_i]._x
        << " " << S._cnState[v_i]._y
        << " " << S._cnState[v_i]._z;
    for (double s : S._cnState[v_i]._s)
    {
      out << " " << s;
    }
    out << std::endl;
  }
  
  // output labels
  out << S.label(0);
  for (int i = 1; i < S.k(); ++i)
  {
    out << " " << S.label(i);
  }
  out << std::endl;

  return out;
}
  
std::istream& operator>>(std::istream& in, StateTree& S)
{
  typedef std::vector<std::string> StringVector;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of(" \t"));
  
  IntVector pi(s.size(), -1);
  for (int i = 0; i < s.size(); ++i)
  {
    ss >> pi[i];
  }
  
  S = StateTree(pi);
  
  for (int i = 0; i < S.k(); ++i)
  {
    getline(in, line);
    ss.clear();
    ss.str(line);
    
    Node v_i = S.node(i);
    StateTree::CopyNumberState& cnState_i = S._cnState[v_i];
    
    boost::split(s, line, boost::is_any_of(" "));
    
    cnState_i._x = boost::lexical_cast<int>(s[0]);
    cnState_i._y = boost::lexical_cast<int>(s[1]);
    cnState_i._z = boost::lexical_cast<int>(s[2]);

    for (int j = 3; j < s.size(); ++j)
    {
      cnState_i._s.push_back(boost::lexical_cast<double>(s[j]));
    }
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  for (int i = 0; i < S.k(); ++i)
  {
    ss >> S._label[S._stateToNode[i]];
  }
  
  return in;
}

std::ostream& operator<<(std::ostream& out, const StateTreeVector& vecS)
{
  out << vecS.size() << " #state trees" << std::endl;
  for (const StateTree& S : vecS)
  {
    out << S;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, StateTreeVector& vecS)
{
  g_lineNumber = 0;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  int nrStateTrees = -1;
  ss >> nrStateTrees;
  
  if (nrStateTrees < 1)
  {
    throw std::runtime_error(getLineNumber() + "Error: incorrect number of state trees");
  }
  
  for (int i = 0; i < nrStateTrees; ++i)
  {
    StateTree S;
    in >> S;
    vecS.push_back(S);
  }
  
  return in;
}
