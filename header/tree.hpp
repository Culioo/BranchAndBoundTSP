//
// Created by adyck on 1/15/18.
//

#ifndef BRANCHANDBOUNDTSP_TREE_HPP
#define BRANCHANDBOUNDTSP_TREE_HPP

#include <cstdlib>
#include <vector>
#include <stdexcept>
#include "graph.hpp"

namespace TSP {
using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;


class Tree : public TSP::Graph {
 public:
  Tree(NodeId root, NodeId size) : TSP::Graph(size) {
    _nodes.push_back(root);
    _root = root;
  }

  TSP::NodeId predecessor(TSP::NodeId id) const {
    auto k = this->node(id).neighbors().at(0);
    if (k == TSP::invalid_node_id)
      throw std::runtime_error("No neighbors on this node!");
    else
      return k;
  }

  const std::vector<NodeId> &nodes() const {
    return _nodes;
  }

  const NodeId &root() const {
    return _root;
  }
 private:
  NodeId _root;
  std::vector<EdgeId> _edges;
  std::vector<NodeId> _nodes;
  size_type num_edges;
  size_type num_nodes;
};

class OneTree {
  OneTree(NodeId root, size_t size) : MST(0, size - 1), _root(root) {
    _nodes.push_back(root);
    for (NodeId node = 0; node )
  }

 private:
  Tree MST;
  NodeId _root;
  std::vector<EdgeId> _edges;
  std::vector<NodeId> _nodes;
  size_type num_edges;
  size_type num_nodes;
};

};


#endif //BRANCHANDBOUNDTSP_TREE_HPP
