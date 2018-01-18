//
// Created by adyck on 1/15/18.
//

#ifndef BRANCHANDBOUNDTSP_TREE_HPP
#define BRANCHANDBOUNDTSP_TREE_HPP

#include<cstdlib>
#include <vector>
#include <stdexcept>
#include "graph.hpp"
#include "util.hpp"

namespace TSP {
using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;

//class Tree : public TSP::Graph {
// public:
//  Tree(NodeId root, NodeId size) : TSP::Graph(size) {
//      _nodes.push_back(root);
//      _root = root;
//  }
//
//  TSP::NodeId predecessor(TSP::NodeId id) const {
//      auto k = this->node(id).neighbors().at(0);
//      if (k == TSP::invalid_node_id)
//          throw std::runtime_error("No neighbors on this node!");
//      else
//          return k;
//  }
//
//  const std::vector<NodeId> &nodes() const {
//      return _nodes;
//  }
//
//  const NodeId &root() const {
//      return _root;
//  }
// private:
//  NodeId _root;
//  std::vector<NodeId> _nodes;
//};

class OneTree {
 public:
  OneTree(size_t size) : size(size), //MST(0, size - 1)
                                      _edges(0) {
      _nodes.push_back(Node()); //root = 0 !!
      for (NodeId node = 1; node < size; node++)
          _nodes.push_back(Node());
      num_nodes = _nodes.size();
      num_edges = 0;
  }

//  OneTree(const OneTree &obj) {
//      _root = obj._root;
//      size = obj.size;
//
//      _nodes = obj._nodes;
//      _edges = obj._edges;
//
//      num_edges = obj.num_edges;
//      num_nodes = obj.num_nodes;
//  }

  void clear() {
      _edges.clear();
      num_edges = 0;
      for (auto &node : _nodes) {
          node._neighbors.clear();
      }
      _root._neighbors.clear();
  }
  void add_edge(NodeId i, NodeId j) {
      if (i >= num_nodes || j >= num_nodes)
          throw std::runtime_error("Index out of range while adding edge to 1-tree");

      _edges.push_back(to_EdgeId(i, j, size));
      _nodes.at(i).add_neighbor(j);
      _nodes.at(j).add_neighbor(i);
      num_edges++;
  }

  const size_type &get_num_edges() const {
      return num_edges;
  }
  const std::vector<EdgeId> &get_edges() const {
      return _edges;
  }

  const std::vector<Node> &get_nodes() const {
      return _nodes;
  }

  const Node &get_node(NodeId id) const {
      return _nodes.at(id);
  }

 private:
//  Tree MST;
  Node _root;
  size_type size;

  std::vector<Node> _nodes;
  std::vector<EdgeId> _edges;

  size_type num_edges;
  size_type num_nodes;
};

};

#endif //BRANCHANDBOUNDTSP_TREE_HPP
