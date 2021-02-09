/*
 *  statetree.h
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef STATETREE_H
#define STATETREE_H

#include "utils.h"
#include "readmatrix.h"

// Forward declaration
class StateTree;

/// Vector of state trees
typedef std::vector<StateTree> StateTreeVector;

/// This class models a state tree
class StateTree
{
private:
  /// This struct models an SNV-specific copy number state
  struct CopyNumberState
  {
    /// Default constructor
    CopyNumberState()
      : _x(-1)
      , _y(-1)
      , _z(-1)
      , _s()
    {
    }
    
    /// Number of maternal copies
    int _x;
    /// Number of paternal copies
    int _y;
    /// Number of mutated copies
    int _z;
    /// SNV-specific copy-number mixture proportions
    DoubleVector _s;
  };
  
public:
  /// Default constructor
  StateTree();
  
  /// Constructor
  ///
  /// @param k Number of states
  StateTree(int k);
  
  /// Copy constructor
  ///
  /// @param other Other state tree
  StateTree(const StateTree& other);
  
  /// Assignment operator
  ///
  /// @param other Other state tree
  StateTree& operator=(const StateTree& other);
  
  /// Constructor
  ///
  /// @param pi Vector of parental states
  StateTree(const IntVector& pi);
  
  /// Write summarize
  ///
  /// @param R Read matrix
  /// @param sol State tree vector
  /// @param out Output stream
  static void writeSummary(const ReadMatrix& R,
                           const StateTreeVector& sol,
                           std::ostream& out);
 
  /// Return number of distinct states
  int k() const
  {
    return _k;
  }
  
  int getNrSamples() const
  {
    return _cnState[_root]._s.size();
  }
  
  /// Get mixture proportions
  ///
  /// @param i State
  DoubleVector getMixtureProportion(int i) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    
    return _cnState[v_i]._s;
  }
  
  void setMixtureProportion(int i, double s)
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    
    _cnState[v_i]._s.push_back(s);
  }
  
  void incrementMixtureProportion(int p, int i, double s)
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    
    _cnState[v_i]._s[p] += s;
  }
  
  void resetMixtureProportions(int nrSamples)
  {
    for (int j = 0; j < _k; ++j)
    {
      Node v_j = _stateToNode[j];
      if (v_j != lemon::INVALID)
      {
        _cnState[v_j]._s = DoubleVector(nrSamples, 0);
      }
    }
  }
  
  void resetMixtureProportions(int p, double dcf)
  {
    const Node v_i = mutationNode(_root);
    const int i = _nodeToState[v_i];
    const Node v_pi_i = _stateToNode[parent(i)];
    
    const double marginal_sum = _cnState[v_i]._s[p] + _cnState[v_pi_i]._s[p];
    
    double sum = 0;
    for (int j = 0; j < _k; ++j)
    {
      if (i == j) continue;
      Node v_j = _stateToNode[j];
      if (v_j == lemon::INVALID)
        continue;
      
      if (isAncestor(i, j))
      {
        sum += _cnState[v_j]._s[p];
      }
    }
    
    _cnState[v_i]._s[p] = dcf - sum;
    _cnState[v_pi_i]._s[p] = marginal_sum - _cnState[v_i]._s[p];
  }
  
  /// Set copy number state
  ///
  /// @param i State
  /// @param x Number of maternal copies
  /// @param y Number of paternal copies
  /// @param z Number of mutated copies
  void setCnState(int i, int x, int y, int z)
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    _cnState[v_i]._x = x;
    _cnState[v_i]._y = y;
    _cnState[v_i]._z = z;
  }
  
  /// Get copy number state
  ///
  /// @param i State
  /// @param x Number of maternal copies
  /// @param y Number of paternal copies
  /// @param z Number of mutated copies
  void getCnState(int i, int& x, int& y, int& z) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    x = _cnState[v_i]._x;
    y = _cnState[v_i]._y;
    z = _cnState[v_i]._z;
  }
  
  /// Return parent state; -1 is returned if i is the root state,
  /// -2 is returned if state is absent from the state tree
  ///
  /// @param i State
  int parent(int i) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    Arc a = InArcIt(_S, v_i);
    if (a == lemon::INVALID)
    {
      if (v_i == _root)
      {
        return -1;
      }
      else
      {
        return -2;
      }
    }
    else
    {
      return _nodeToState[_S.source(a)];
    }
  }
  
  /// Return number of vertices
  int numVertices() const
  {
    return lemon::countNodes(_S);
  }
  
  /// Return whether specified state is present
  ///
  /// @param i State
  bool isPresent(int i) const
  {
    assert(0 <= i && i < _k);
    return parent(i) != -2;
  }
  
  /// Return whether specified states are in a parent-child relationship
  ///
  /// @param i Parent state
  /// @param j Child state
  bool isParent(int i, int j) const
  {
    assert(0 <= i && i < _k);
    assert(1 <= j && j < _k);

    Node v_i = _stateToNode[i];
    Node v_j = _stateToNode[j];
    
    Arc a = InArcIt(_S, v_j);
    if (a == lemon::INVALID)
    {
      return false;
    }
    else
    {
      return v_i == _S.source(a);
    }
  }
  
  /// Return whether specified states are siblings
  ///
  /// @param i State
  /// @param j State
  bool areSibblings(int i, int j) const
  {
    assert(1 <= i && i < _k);
    assert(1 <= j && j < _k);
    
    Node v_i = _stateToNode[i];
    Node v_j = _stateToNode[j];
    
    Arc a1 = InArcIt(_S, v_i);
    if (a1 == lemon::INVALID)
      return false;
    
    Node v_pi_i = _S.source(a1);
    
    Arc a2 = InArcIt(_S, v_j);
    if (a2 == lemon::INVALID)
      return false;
    Node v_pi_j = _S.source(a2);
    
    return v_pi_i == v_pi_j;
  }
  
  /// Return whether specified states are in an ancestor-descendant relationship
  ///
  /// @param i Ancestor state
  /// @param j Descendant state
  bool isAncestor(int i, int j) const
  {
    assert(0 <= i && i < _k);
    assert(0 <= j && j < _k);
    
    Node v_i = _stateToNode[i];
    
    Node v = _stateToNode[j];
    while (v != v_i && v != _root)
    {
      v = _S.source(InArcIt(_S, v));
    }
    
    return v == v_i;
  }
  
  /// Return whether specified states are in a descendant-ancestor relationship
  ///
  /// @param i Descendant state
  /// @param j Ancestor state
  bool isDescendant(int i, int j) const
  {
    return isAncestor(j, i);
  }
  
  /// Return whether specified states occur in distinct branches
  ///
  /// @param i State
  /// @param j State
  bool isIncomparable(int i, int j) const
  {
    return !isAncestor(i, j) && !isAncestor(j, i);
  }
  
  /// Return set of descendant states of specified state
  ///
  /// @param i State
  const IntSet& D(int i) const
  {
    assert(0 <= i && i < _k);
    return _D[_stateToNode[i]];
  }
  
  /// Return label of specified state
  ///
  /// @param i State
  const std::string& label(int i) const
  {
    assert(0 <= i && i < _k);
    return _label[_stateToNode[i]];
  }
  
  /// Set label of specified state
  ///
  /// @param i State
  /// @param l Label
  void setLabel(int i, const std::string& l)
  {
    _label[_stateToNode[i]] = l;
  }
  
  /// Return state tree
  const Digraph& S() const
  {
    return _S;
  }
  
  /// Return the state of the specified node
  ///
  /// @param v_i Node
  int state(Node v_i) const
  {
    assert(v_i != lemon::INVALID);
    return _nodeToState[v_i];
  }
  
  /// Return the node of the specified state
  ///
  /// @param i State
  Node node(int i) const
  {
    assert(0 <= i && i < _k);
    return _stateToNode[i];
  }
  
  /// Return root state
  int rootState() const
  {
    return _nodeToState[_root];
  }
  
  /// Write edge list
  ///
  /// @param out Output stream
  void writeEdgeList(std::ostream& out) const;
  
  /// Write state tree in DOT format
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  void writeClusterDOT(double gamma,
                       int t,
                       std::ostream& out) const;
  
  /// Write state tree in DOT format
  ///
  /// @param m Number of samples
  /// @param out Output stream
  void writeDOT(const int m,
                std::ostream& out) const;
  
  /// Return number of maternal copies of specified state
  ///
  /// @param i State
  int x(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._x;
  }
  
  /// Return number of paternal copies of specified state
  ///
  /// @param i State
  int y(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._y;
  }

  /// Return number of mutated copies of specified state
  ///
  /// @param i State
  int z(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._z;
  }
  
  /// Return SNV-specific copy number mixture proportion
  ///
  /// @param i State
  /// @param p Sample
  double s(int i, int p) const
  {
    Node v_i = node(i);
    assert(0 <= p && p < _cnState[v_i]._s.size());
    return _cnState[v_i]._s[p];
  }
  
  /// Return the unique state whose incoming edge is a mutation edge
  int mutationState() const
  {
    Node mut_node = mutationNode(_root);
    if (mut_node == lemon::INVALID)
    {
      return -1;
    }
    else
    {
      return _nodeToState[mut_node];
    }
  }
  
  /// Return variant alllele frequency of the specified sample
  ///
  /// @param p Sample
  double vaf(int p) const
  {
    double numerator = 0;
    double denominator = 0;
    for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
    {
      const int i = _nodeToState[v_i];
      if (isPresent(i))
      {
        numerator += _cnState[v_i]._z * _cnState[v_i]._s[p];
        denominator += (_cnState[v_i]._x + _cnState[v_i]._y) * _cnState[v_i]._s[p];
      }
    }
    return numerator / denominator;
  }
  
  /// Return descendant cell fraction of the specified sample
  ///
  /// @param p Sample
  double dcf(int p) const
  {
    double res = 0;
    
    int i = mutationState();
    assert(i != -1);
    
    for (int j : D(i))
    {
      Node v_j = _stateToNode[j];
      res += _cnState[v_j]._s[p];
    }
    
    return res;
  }
  
  /// Return cell fraction of the specified sample
  ///
  /// @param p Sample
  double cf(int p) const
  {
    double res = 0;
    
    for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
    {
      const int i = _nodeToState[v_i];
      if (isPresent(i) && _cnState[v_i]._z > 0)
      {
        res += _cnState[v_i]._s[p];
      }
    }

    return res;
  }
  
  /// Return maximum likelihood CCF given vaf
  ///
  /// @param p Sample
  /// @param vaf
  double maxLikelihoodCCF(int p,
                          double vaf) const;
  
  /// Return maximum likelihood DCF given vaf
  ///
  /// @param p Sample
  /// @param vaf
  double maxLikelihoodDCF(int p,
                          double vaf) const;
  
  /// Generate mixture proportions for specified sample and dcf
  ///
  /// @param f DCF
  /// @param purity Purity
  void sampleMixtureProportions(double f, double purity);
  
  /// Clear mixture proportions
  void clearMixtureProportions()
  {
    for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
    {
      _cnState[v_i]._s.clear();
    }
  }
  
  /// Check whether mixture proportions are ok
  ///
  /// @param m Number of samples
  bool okMixtureProportions(int m) const
  {
    for (NodeIt v_i(_S); v_i != lemon::INVALID; ++v_i)
    {
      const int i = _nodeToState[v_i];
      const CopyNumberState& cnState_i = _cnState[v_i];
      if (isPresent(i))
      {
        if (cnState_i._s.size() != m)
        {
          return false;
        }
        
        bool ok = false;
        for (int p = 0; p < m; ++p)
        {
          ok |= g_tol.nonZero(cnState_i._s[p]);
        }
        if (!ok)
        {
          return false;
        }
      }
    }
    
    return true;
  }

private:
  typedef Digraph::NodeMap<CopyNumberState> CopyNumberStateNodeMap;
  
  Node computeMaxLikelihoodStateProportions(DoubleNodeMap& s,
                                            int p,
                                            double vaf) const;
  
private:
  /// Number of states
  int _k;
  /// State tree
  Digraph _S;
  /// Root node
  Node _root;
  /// node labeling
  StringNodeMap _label;
  /// State to node mapping
  NodeVector _stateToNode;
  /// node to state mapping
  IntNodeMap _nodeToState;
  /// Descendant state set node map
  IntSetNodeMap _D;
  /// Copy number state node map
  CopyNumberStateNodeMap _cnState;

  /// Initialize descendant sets of subtree rooted at the specified node
  ///
  /// @param v_i node
  void initD(Node v_i);
  
  /// Initialize state tree using specificied parental state vector
  ///
  /// @param pi Vector of parental states
  void init(const IntVector& pi);
  
  /// Identify mutation node recursively from v_i
  ///
  /// @param v_i Node
  Node mutationNode(Node v_i) const
  {
    if (_cnState[v_i]._z > 0)
    {
      return v_i;
    }
    else
    {
      for (OutArcIt a_ij(_S, v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = _S.target(a_ij);
        Node mut_node = mutationNode(v_j);
        if (mut_node != lemon::INVALID)
        {
          return mut_node;
        }
      }
      
      return lemon::INVALID;
    }
  }
  
  friend std::ostream& operator<<(std::ostream& out, const StateTree& S);
  friend std::istream& operator>>(std::istream& in, StateTree& S);
  friend std::ostream& operator<<(std::ostream& out, const StateTreeVector& vecS);
  friend std::istream& operator>>(std::istream& in, StateTreeVector& vecS);
};

/// Output state tree
///
/// @param out Output stream
/// @param S State tree
std::ostream& operator<<(std::ostream& out, const StateTree& S);

/// Input state tree
///
/// @param in Input stream
/// @param S State tree
std::istream& operator>>(std::istream& in, StateTree& S);

/// Output state tree vector
///
/// @param out Output stream
/// @param vecS State tree vector
std::ostream& operator<<(std::ostream& out, const StateTreeVector& vecS);

/// Input state tree vector
///
/// @param in Input stream
/// @param vecS State tree vector
std::istream& operator>>(std::istream& in, StateTreeVector& vecS);

#endif // STATETREE_H
