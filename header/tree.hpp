//
// Created by adyck on 1/15/18.
//

#ifndef BRANCHANDBOUNDTSP_TREE_HPP
#define BRANCHANDBOUNDTSP_TREE_HPP

#include <cstdlib>
#include <vector>
#include "graph.hpp"

namespace TSP {
using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;
//class Tree {
// public:
//  Tree () {
//
//  }
//
//  Tree(NodeId root) : _root(root) {
//
//  }
//
//  void add
//
// private:
//  NodeId  _root;
//  std::vector<NodeId> _nodes;
//  std::vector<EdgeId> _edges;


class Tree : TSP::Graph {
 public:
  Tree() {

  }

  void add_edge(NodeId i, NodeId j) {

  }

 private:
  std::vector<EdgeId> _edges;
  std::vector<NodeId> _nodes;
  size_type num_edges;
  size_type num_nodes;
};



};


#endif //BRANCHANDBOUNDTSP_TREE_HPP
