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
#include "util.hpp"


// TODO: Split header and cpp methods
namespace TSP {

using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;

EdgeId to_EdgeId(NodeId i, NodeId j, size_type N) {
    if (i == j)
        throw std::runtime_error("Loops are not contained in this instance");
    if (i > j)
        std::swap(i, j);
    return (i) * (N) + j;
}

void to_NodeId(EdgeId e, NodeId &i, NodeId &j, size_type N) {
    j = e % (N);
    i = (e - j) / (N);
}

template<class coord_type, class dist_type>
class BranchingNode;

template<class coord_type, class dist_type>
class Instance;

template<class coord_type, class dist_type>
void compute_minimal_1_tree(std::vector<EdgeId> &tree,
                            std::vector<NodeId> &tree_deg,
                            const std::vector<dist_type> lambda,
                            const Instance<coord_type, dist_type> &tsp,
                            const BranchingNode<coord_type, dist_type> &bn) {
    // Initializing union find structure, i and j nodeIDs and a sorted vector
    // which will contain the edges in an sorted order
    Union_Find uf(tsp.size());
    NodeId i = 0, j = 0;
    std::vector<EdgeId> sorted(tsp.num_edges());
    for (size_type index = 0; index < tsp.num_edges(); index++)
        sorted[index] = index;

    for (auto const &el : bn.get_required()) {
        tree.push_back(el);
        to_NodeId(el, i, j, tsp.size());
        uf._union(i, j);
        tree_deg.at(i)++, tree_deg.at(j)++;
    }

    std::vector<dist_type> mod_weights;
    mod_weights.reserve(tsp.num_edges());
    for (size_t lauf = 0; lauf < tsp.weights().size(); lauf++) {
        NodeId index1 = 0, index2 = 0;
        to_NodeId(lauf, index1, index2, tsp.size());
        mod_weights.push_back(tsp.weight(lauf) + lambda[index1] + lambda[index2]);
        if(bn.check_required(lauf)){
            *mod_weights.rbegin() = -std::numeric_limits<dist_type>::max();
        } else
        if (bn.check_forbidden(lauf)) {
            *mod_weights.rbegin() = std::numeric_limits<dist_type>::max();
        }
    }
    std::stable_sort(sorted.begin(), sorted.end(),
                     [&](EdgeId ind1, EdgeId ind2) {
                       return (mod_weights[ind1] <
                           mod_weights[ind2]);
                     });
    std::vector<EdgeId> candidates;
    candidates.reserve(tsp.size());
    for (auto it = sorted.begin() ; (it != sorted.end()) ; it++) {
        to_NodeId(*it, i, j, tsp.size());
        if ((i == 0) && (0 != j) && ( candidates.size() < 2 ))
            candidates.push_back(*it);
        if ((i != j) && (i != 0) && (0 != j) && (tree.size() != tsp.size() - 1)) {
            if (uf._find(i) != uf._find(j)) {
                tree.push_back(*it);
                uf._union(i, j);
                tree_deg.at(i)++, tree_deg.at(j)++;
            }
        }
    }
    tree.push_back(candidates.at(0));
    to_NodeId(candidates.at(0), i, j, tsp.size());
    tree_deg.at(i)++, tree_deg.at(j)++;

    tree.push_back(candidates.at(1));
    to_NodeId(candidates.at(1), i, j, tsp.size());
    tree_deg.at(i)++, tree_deg.at(j)++;
    return;
}

template<class coord_type, class dist_type>
dist_type Held_Karp(const TSP::Instance<coord_type, dist_type> &tsp,
                    std::vector<dist_type> lambda,
                    std::vector<EdgeId> &tree,
                    std::vector<NodeId> &tree_deg,
                    const BranchingNode<coord_type, dist_type> &bn) {
    size_type n = tsp.size();
    std::vector<dist_type> sol_vector;
    std::vector<dist_type> lambda_max(lambda), lambda_tmp(lambda);
    size_t max_el = 0;
    std::vector<EdgeId> tree_max(tsp.size()), tree_deg_max(tree_deg);
    for (size_t i = 0; i < std::ceil(n * n / 50) + n + 15; i++) {
        tree.clear();
        for (auto it = tree_deg.begin(); it != tree_deg.end(); it++)
            *it = 0;
        compute_minimal_1_tree<coord_type, dist_type>(tree, tree_deg, lambda_tmp, tsp, bn);
        if (i == 0) {
            tree_max = tree;
            tree_deg_max = tree_deg;
        }

        dist_type sum = 0, sum2 = 0;
        for (size_t k1 = 0; k1 < tree.size(); k1++)
            sum += tsp.weight(tree[k1]);
        for (size_t k2 = 0; k2 < tree_deg.size(); k2++)
            sum2 += (tree_deg[k2] - 2.) * lambda_tmp[k2];
        sol_vector.push_back(sum + sum2);

        if (i > 0)
            if (sol_vector[max_el] < sol_vector[i]) {
                lambda_max = lambda_tmp;
                max_el = i;
                tree_max = tree;
                tree_deg_max = tree_deg;
            }

        for (size_t j = 0; j < lambda_tmp.size(); j++) {
            lambda_tmp[j] += 1 * (tree_deg[j] - 2.);  // TODO: replace 1 by t_i
        }
    }
    lambda = lambda_max;
    tree = tree_max;
    tree_deg = tree_deg_max;
    return *std::max_element(sol_vector.begin(), sol_vector.end());
}

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
      return _weights[id];
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
  BranchingNode(const std::vector<EdgeId> &req,
                const std::vector<EdgeId> &forb,
                const Instance<coord_type, dist_type> &tsp,
                const std::vector<dist_type> &lambda
  ) : required(req), forbidden(forb), lambda(lambda), tree(tsp.size()), tree_deg(tsp.size()) {
      this->HK = Held_Karp(tsp, this->lambda, this->tree, this->tree_deg, *this);
  }

  BranchingNode(const std::vector<EdgeId> &req,
                const std::vector<EdgeId> &forb,
                const Instance<coord_type, dist_type> &tsp,
                const std::vector<dist_type> &lambda,
                EdgeId e1
  ) : required(req), forbidden(forb), lambda(lambda), tree(tsp.size()), tree_deg(tsp.size()) {
      forbidden.push_back(e1);
      this->HK = Held_Karp(tsp, this->lambda, this->tree, this->tree_deg, *this);

  }

  BranchingNode(const std::vector<EdgeId> &req,
                const std::vector<EdgeId> &forb,
                const Instance<coord_type, dist_type> &tsp,
                const std::vector<dist_type> &lambda,
                EdgeId e1,
                EdgeId e2
  ) : required(req), forbidden(forb), lambda(lambda), tree(tsp.size()), tree_deg(tsp.size()) {
      this->HK = Held_Karp(tsp, this->lambda, this->tree, this->tree_deg, *this);
//      required.push_back(e1);
      forbidden.push_back(e2);

  }

  BranchingNode(const std::vector<EdgeId> &req,
                const std::vector<EdgeId> &forb,
                const Instance<coord_type, dist_type> &tsp,
                const std::vector<dist_type> &lambda,
                EdgeId e1,
                EdgeId e2,
                bool both
  ) : required(req), forbidden(forb), lambda(lambda), tree(tsp.size()), tree_deg(tsp.size()) {
      this->HK = Held_Karp(tsp, this->lambda, this->tree, this->tree_deg, *this);
      if (both) {
          required.push_back(e1);
          required.push_back(e2);
      }

  }

  bool check_required(EdgeId id) const {
      return (std::find(required.begin(), required.end(), id) != required.end());
  }

  bool check_forbidden(EdgeId id) const {
      return (std::find(forbidden.begin(), forbidden.end(), id) != forbidden.end());
  }

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
  const std::vector<EdgeId> &get_tree() const {
      return tree;
  }
  std::vector<EdgeId> &get_tree() {
      return tree;
  }

  const std::vector<EdgeId> &get_tree_deg() const {
      return tree;
  }
  std::vector<EdgeId> &get_tree_deg() {
      return tree;
  }

  void set_HK(dist_type c) {
      this->HK = c;
  }
  dist_type get_HK() {
      return this->HK;
  }

  bool tworegular() {
      for (auto it = tree_deg.begin(); it != tree_deg.end(); it++)
          if (*it != 2)
              return false;
      return true;
  }

 private:
  friend Instance<coord_type, dist_type>;
  std::vector<EdgeId> required;
  std::vector<EdgeId> forbidden;
  std::vector<dist_type> lambda;
  std::vector<EdgeId> tree;
  std::vector<NodeId> tree_deg;
  dist_type HK;
};
}

#include "tsp_impl.hpp"

#endif // BRANCHANDBOUNDTSP_TSP_HPP