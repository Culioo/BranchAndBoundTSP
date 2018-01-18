// modified by Alex Dyck and Leon Sievers 10/31/17.
#ifndef GRAPH_HPP
#define GRAPH_HPP

/**
   @file graph.hpp

   @brief This file provides a simple class @c Graph to model unweighted undirected graphs.
**/

#include <cstddef> // std::size_t
#include <iosfwd> // std::ostream fwd declare
#include <limits>
#include <vector>

namespace TSP // for Edmonds
{

using size_type = std::size_t;

//! Always use these typedefs to identify nodes by their numbers.
//! Use different typedefs for internal standard indexing starting with 0
//! and DIMACS-based indexing starting with 1.
//! Even better, one could use strong typedefs (i.e. use an enum class
//! or a custom struct providing cast operator to \c size_type s.t.
//! mixup of index types is avoided.
//! This is not done here for simplicity.
using NodeId = size_type;
using TSPlibId = size_type;

/** Useful constant different from the id of any actual node: **/
NodeId constexpr invalid_node_id = std::numeric_limits<NodeId>::max();
TSPlibId constexpr invalid_dimacs_id = std::numeric_limits<TSPlibId>::max();

/**
   Nodes in DIMACS files are counted from 1, but here we count them from 0 so they match their std::vector indices.
   These two trivial functions should help make the transition between the two models clear (instead of just having
   some unexplained -1's and +1's in the middle of the code.
**/
NodeId from_tsplib_id(TSPlibId const tsplib_id); //!< Subtracts 1 (throws if @c dimacs_id is 0)
TSPlibId to_tsplib_id(NodeId const node_id);     //!< Adds 1 (throws if overflow would occur)

/**
   @class Node

   @brief A @c Node stores an array of neighbors (via their ids).

   @note The neighbors are not necessarily ordered, so searching for a specific neighbor takes O(degree)-time.
**/
class Node {
 public:
  typedef std::size_t size_type;

  /** @brief Create an isolated node (you can add neighbors later). **/
  Node() {};   //= default;

  /** @return The number of neighbors of this node. **/
  size_type degree() const;

  /** @return The array of ids of the neighbors of this node. **/
  std::vector<NodeId> const &neighbors() const;
  /**
     @brief Adds @c id to the list of neighbors of this node.
     @warning Does not check whether @c id is already in the list of neighbors (a repeated neighbor is legal, and
     models parallel edges).
     @warning Does not check whether @c id is the identity of the node itself (which would create a loop!).
  **/
  void add_neighbor(NodeId const id);
 private:
  friend class Graph;
  friend class OneTree;
  friend class BranchingTree;

  std::vector<NodeId> _neighbors;
}; // class Node

//BEGIN: Inline section

inline
Node::size_type Node::degree() const {
    return neighbors().size();
}

inline
std::vector<NodeId> const &Node::neighbors() const {
    return _neighbors;
}





//END: Inline section



} // namespace ED

#endif /* GRAPH_HPP */
