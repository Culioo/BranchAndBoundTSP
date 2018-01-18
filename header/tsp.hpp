#ifndef BRANCHANDBOUNDTSP_TSP_HPP
#define BRANCHANDBOUNDTSP_TSP_HPP

#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <climits>
#include <vector>
#include <sstream>
#include <queue>
#include <numeric>
#include <utility>
#include <cassert>
#include "util.hpp"
#include "tree.hpp"

#define EPS 10e-7


// TODO: Split header and cpp methods
namespace TSP {

using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;

template<class coord_type, class dist_type>
class BranchingNode;

template<class coord_type, class dist_type>
class Instance;

template<class coord_type, class dist_type>
class Instance {
 public:
  Instance(const std::string &filename);

  dist_type distance(coord_type x1, coord_type y1, coord_type x2, coord_type y2) {
      return std::lround(std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
  }

  void compute_optimal_tour();

  void print_optimal_length();

  void print_optimal_tour(const std::string &filename);

  size_type size() const {
      return dimension;
  }

  size_type num_edges() const {
      return _weights.size();
  }

  dist_type weight(EdgeId id) const {
      return _weights.at(id);
  }

  const std::vector<dist_type> &weights() const {
      return _weights;
  }

 private:
  std::vector<NodeId> _nodes;
  std::vector<dist_type> _weights;
  size_type dimension;
  std::vector<NodeId> _tour;
};

template<class coord_type, class dist_type>
class BranchingNode {
 public:
  BranchingNode(const Instance<coord_type, dist_type> &tsp
  ) : size(tsp.size()), required(), required_neighbors(size), forbidden(),
      forbidden_neighbors(size), lambda(size, 0), tree(size) {
      HK = Held_Karp(tsp, this->lambda, this->tree, *this, true);
  }

  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1
  ) : size(tsp.size()),
      required(BNode.get_required()),
      required_neighbors(BNode.required_neighbors),
      forbidden(BNode.get_forbidden()),
      forbidden_neighbors(BNode.forbidden_neighbors),
      lambda(BNode.get_lambda()),
      tree(tsp.size()) {

      add_forbidden(e1);
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }

  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1,
                EdgeId e2
  ) : size(tsp.size()),
      required(BNode.get_required()),
      required_neighbors(BNode.required_neighbors),
      forbidden(BNode.get_forbidden()),
      forbidden_neighbors(BNode.forbidden_neighbors),
      lambda(BNode.get_lambda()),
      tree(tsp.size()) {
      add_required(e1);
      add_forbidden(e2);
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }

  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1,
                EdgeId e2,
                bool both_req
  ) : size(tsp.size()),
      required(BNode.get_required()),
      required_neighbors(BNode.required_neighbors),
      forbidden(BNode.get_forbidden()),
      forbidden_neighbors(BNode.forbidden_neighbors),
      lambda(BNode.get_lambda()),
      tree(tsp.size()) {

      if (both_req) {
          add_required(e1);
          add_required(e2);
      }
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }

  bool operator<(const BranchingNode<coord_type, dist_type> &rhs) const;

  EdgeId reverse_edge(EdgeId e, size_type n) const {
      NodeId i = 0, j = 0;
      to_NodeId(e, i, j, n);
      return j * n + i;
  }

  bool is_required(EdgeId id) const {
      return (std::find(required.begin(), required.end(), id) != required.end()
          && std::find(required.begin(), required.end(), reverse_edge(id, size)) != required.end());
  }

  bool is_forbidden(EdgeId id) const {
      return (std::find(forbidden.begin(), forbidden.end(), id) != forbidden.end()
          && std::find(forbidden.begin(), forbidden.end(), reverse_edge(id, size)) != forbidden.end());
  }

  void forbid(NodeId idx, EdgeId e1, EdgeId e2) ;

  void admit(NodeId idx) ;

  bool push_required(EdgeId e);

  void add_required(EdgeId e) ;

  bool push_forbidden(EdgeId e);

  void add_forbidden(EdgeId e);

  const std::vector<EdgeId> &get_required() const {
      return required;
  }
  const std::vector<EdgeId> &get_forbidden() const {
      return forbidden;
  }

  const std::vector<dist_type> &get_lambda() const {
      return lambda;
  }
  std::vector<dist_type> &get_lambda() {
      return lambda;
  }

  const OneTree &get_tree() const {
      return tree;
  }
  OneTree &get_tree() {
      return tree;
  }

  void set_HK(dist_type c) {
      this->HK = c;
  }
  const dist_type get_HK() const {
      return this->HK;
  }

  const std::vector<Node> &get_required_neighbors() const {
      return required_neighbors;
  }

  bool tworegular() {
      for (const auto &el : tree.get_nodes())
          if (el.degree() != 2)
              return false;
      return true;
  }

 private:
  size_type size;
  std::vector<EdgeId> required;
  std::vector<Node> required_neighbors;

  std::vector<EdgeId> forbidden;
  std::vector<Node> forbidden_neighbors;

  std::vector<double> lambda;
  OneTree tree;

  dist_type HK;
};
}

#include "tsp_impl.hpp"

#endif // BRANCHANDBOUNDTSP_TSP_HPP