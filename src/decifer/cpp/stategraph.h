/*
 * stategraph.h
 *
 *  Created on: 11-jan-2016
 *      Author: M. El-Kebir
 */

#ifndef STATEGRAPH_H
#define STATEGRAPH_H

#include <lemon/adaptors.h>
#include <vector>
#include "utils.h"

/// This class models a directed graph containing all state trees
class StateGraph
{
public:
  /// Constructor
  ///
  /// @param max_x Maximum number of maternal copies
  /// @param max_y Maximum number of paternal copies
  StateGraph(int max_x, int max_y);
  
  /// Mutation type
  typedef enum { MUTATION, AMPLIFICATION, DELETION } EdgeType;
 
  /// Edge type map
  typedef Digraph::ArcMap<EdgeType> EdgeTypeArcMap;
  
  /// Write DOT visualization of state graph
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Enumerate state trees and return the number of enumerated state trees
  ///
  /// @param L Copy number states that must be present
  int enumerate(const IntPairSet& L);
  
private:
  typedef std::vector<NodeMatrix> Node3Matrix;
  typedef std::vector<Node3Matrix> Node4Matrix;
  
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;
  
  /// Initialize state graph
  void init();
  
  /// Initialize enumeration by adding root node to T and setting the frontier F
  ///
  /// @param subG Subgraph of state graph
  /// @param T Partial state tree
  /// @param F Frontier
  void init(SubDigraph& subG,
            SubDigraph& T,
            ArcList& F);
  
  /// Report enumerated state tree
  ///
  /// @param L
  /// @param LL
  /// @param T Enumerated state tree
  /// @param V
  bool finalize(const IntPairSet& L,
                const BoolMatrix& LL,
                const SubDigraph& T,
                const IntMatrix& V,
                const Arc& mutationEdge);
  
  /// Write partial state tree in DOT format
  ///
  /// @param T Partial state tree
  /// @param out Output stream
  void writeDOT(const SubDigraph& T,
                std::ostream& out) const;
  
  /// Write partial state tree and frontier in DOT format
  ///
  /// @param T Partial state tree
  /// @param F Frontier
  /// @param out Output stream
  void writeDOT(const SubDigraph& T,
                const ArcList& F,
                std::ostream& out) const;
  
  /// Add edge to G
  ///
  /// @param v Source node
  /// @param type Edge type (mutation, amplification, deletion)
  /// @param x Number of maternal copies of target node
  /// @param y Number of paternal copies of target node
  /// @param xbar Number of mutated maternal copies of target node
  /// @param ybar Number of mutated paternal copies of target node
  void addEdge(Node v, EdgeType type, int x, int y, int xbar, int ybar)
  {
    Node w = _toNode[x][y][xbar][ybar];
    Arc vw = _G.addArc(v, w);
    _type[vw] = type;
  }
  
  /// Return whether proposed partial state tree is valid
  ///
  /// @param L
  /// @param LL
  /// @param T Proposed partial state tree
  /// @param V
  /// @param verticesCopyStates Copy number states of vertices
  /// @param leavesCopyStates Copy number states of leaves
  bool isValid(const IntPairSet& L,
               const BoolMatrix& LL,
               const SubDigraph& T,
               const IntMatrix& V,
               const IntPairSet& verticesCopyStates,
               const IntPairSet& leavesCopyStates) const;
  
  /// Return whether partial state tree is valid
  ///
  /// @param T Partial state tree
  bool isValid(const SubDigraph& T) const;
  
  /// Return whether partial state tree remains valid upon edge of specified edge
  ///
  /// @param T Partial state tree
  /// @param a Arc
  bool isValid(const SubDigraph& T,
               Arc a);
  
  /// Grow partial state tree
  ///
  /// @param L Set of required copy number states
  /// @param LL Set of required copy number states in matrix form
  /// @param G Graph
  /// @param T Current partial state tree
  /// @param F Frontier of edges to expand into
  /// @param V Nodes of G in matrix form
  /// @param mutationEdge Mutation edge
  /// @param verticesCopyStates Copy number states of vertices
  /// @param leavesCopyStates Copy number states of leaves
  void grow(const IntPairSet& L,
            const BoolMatrix& LL,
            SubDigraph& G,
            SubDigraph& T,
            ArcList& F,
            IntMatrix& V,
            Arc& mutationEdge,
            IntPairSet& verticesCopyStates,
            IntPairSet& leavesCopyStates);
  
public:
  struct CnaTriple
  {
    CnaTriple()
      : _x(-1)
      , _y(-1)
      , _z(-1)
    {
    }
    
    CnaTriple(int x, int y, int z)
      : _x(x)
      , _y(y)
      , _z(z)
    {
    }
    
    bool operator==(const CnaTriple& other) const
    {
      return _x == other._x && _y == other._y && _z == other._z;
    }
    
    bool operator!=(const CnaTriple& other) const
    {
      return !(_x == other._x && _y == other._y && _z == other._z);
    }
    
    int _x;
    int _y;
    int _z;
  };

  typedef std::set<CnaTriple> CnaTripleSet;
  typedef CnaTripleSet::const_iterator CnaTripleSetIt;
  typedef std::vector<CnaTriple> CnaTripleVector;
  
  typedef std::pair<CnaTriple, CnaTriple> StateEdge;
  typedef std::set<StateEdge> StateEdgeSet;
  typedef StateEdgeSet::const_iterator StateEdgeSetIt;
  
  typedef std::set<StateEdgeSet> StateEdgeSetSet;
  typedef StateEdgeSetSet::const_iterator StateEdgeSetSetIt;
  typedef StateEdgeSetSet::iterator StateEdgeSetSetNonConstIt;

private:
  typedef Digraph::NodeMap<CnaTriple> CnaTripleNodeMap;
  
private:
  const int _maxCopyNumberX;
  const int _maxCopyNumberY;
  
  Digraph _G;
  IntNodeMap _x;
  IntNodeMap _y;
  IntNodeMap _xbar;
  IntNodeMap _ybar;
  EdgeTypeArcMap _type;
  
  Node4Matrix _toNode;
  Node _root;
  
  StateEdgeSetSet _result;
  
public:
  static const StateEdgeSetSet& getStateTrees(const IntPairSet& L,
                                              int xy_max_c);
  
  static void writeStateTrees(std::ostream& out);
  
  static void readStateTrees(std::istream& in);
  
  static void partition(const StateEdgeSet& S,
                        const CnaTriple& vertex,
                        bool& preMutationFlag,
                        CnaTriple& mutation,
                        CnaTripleSet& preMutation,
                        CnaTripleSet& postMutation);
  
  static void partition(const StateEdgeSet& S,
                        CnaTripleSet& preMutation,
                        CnaTripleSet& postMutation,
                        CnaTriple& mutation);

  static void print(const StateEdgeSet& S,
                    std::ostream& out)
  {
    for (StateEdgeSetIt it = S.begin(); it != S.end(); ++it)
    {
      out << "( (" << it->first._x << "," << it->first._y << "," << it->first._z << ") , ("
          << it->second._x << "," << it->second._y << "," << it->second._z << ") )  ";
    }
  }
  
  static void print(const StateEdgeSetSet& setS,
                    std::ostream& out)
  {
    for (StateEdgeSetSetIt it1 = setS.begin(); it1 != setS.end(); ++it1)
    {
      const StateEdgeSet& S = *it1;
      print(S, std::cout);
    }
  }
  
public:
  typedef std::map<IntPairSet, StateEdgeSetSet> Dictionary;
  
private:
  static Dictionary _dict;
};
  
bool operator<(const StateGraph::CnaTriple& lhs, const StateGraph::CnaTriple& rhs);

std::ostream& operator<<(std::ostream& out,
                         const StateGraph::Dictionary& dict);

std::istream& operator>>(std::istream& in,
                         StateGraph::Dictionary& dict);

#endif // STATEGRAPH_H
