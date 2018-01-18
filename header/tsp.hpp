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
void compute_minimal_1_tree(OneTree &tree,
                            const std::vector<double> &lambda,
                            const Instance<coord_type, dist_type> &tsp,
                            const BranchingNode<coord_type, dist_type> &BNode) {
    /*
    // Initializing union find structure, i and j nodeIDs and a sorted vector
    // which will contain the edges in an sorted order
    Union_Find uf(tsp.size());
    NodeId i = 0, j = 0;

    std::vector<EdgeId> sorted(tsp.num_edges());
    for (size_type index = 0; index < tsp.num_edges(); index++)
        sorted.at(index) = index;

    //compute modified weights c_\lambda and set the weight of required edges
    // to -inf and for forbidden edges to +inf
    std::vector<dist_type> mod_weights(tsp.num_edges(), 0);
    for (const auto &el : BNode.get_forbidden())
        mod_weights.at(el) = std::numeric_limits<dist_type>::max();
    for (const auto &el : BNode.get_required())
        mod_weights.at(el) = std::numeric_limits<dist_type>::min();


    for (size_t edge = 0; edge < tsp.weights().size(); edge++) {
        if (mod_weights.at(edge) == 0) {//double comparison with == is bad, but let's make it a TODO
            NodeId v = 0, w = 0;
            to_NodeId(edge, v, w, tsp.size());
            mod_weights.at(edge) = tsp.weight(edge) + lambda[v] + lambda[w];
        }
    }
    //Sort the array to compute an 1-tree via iterating over the edges
    std::stable_sort(sorted.begin(), sorted.end(),
                     [&](EdgeId ind1, EdgeId ind2) {
                       return (mod_weights.at(ind1) <
                           mod_weights.at(ind2));
                     });

    std::vector<NodeId> candidates;
    candidates.reserve(tsp.size());
    for (auto it = sorted.begin(); it != sorted.end(); it++) {
        to_NodeId(*it, i, j, tsp.size());
        if ((i == 0) && (0 != j) && (candidates.size() < 2))
            candidates.push_back(j);
        if ((i != j) && (i != 0) && (0 != j) && (tree.get_num_edges() < tsp.size()-1 )) {
            EdgeId edge = to_EdgeId(i,j,tsp.size());
            if ((uf._find(i) != uf._find(j))) {

                assert(!BNode.is_forbidden(edge));
                tree.add_edge(i, j);
                uf._union(i, j);
            }
        }
        if(candidates.size() == 2 && tree.get_num_edges() == tsp.size() - 2) break;
    }
    assert(tree.get_num_edges() == tsp.size() - 2);
    // This should be okay since we forbid all edges where deg >= 3 could happen
    assert(!BNode.is_forbidden(to_EdgeId(0,candidates.at(0),tsp.size())));
    assert(!BNode.is_forbidden(to_EdgeId(0,candidates.at(1),tsp.size())));
    tree.add_edge(0, candidates.at(0));
    tree.add_edge(0, candidates.at(1));
     */

    //compute modified weights c_\lambda and set the weight of required edges
    // to -inf and for forbidden edges to +inf
    std::vector<dist_type> mod_weights(tsp.num_edges(), 0);


    for (size_t edge = 0; edge < tsp.weights().size(); edge++) {
            NodeId v = 0, w = 0;
            to_NodeId(edge, v, w, tsp.size());
            mod_weights.at(edge) = tsp.weight(edge) + lambda[v] + lambda[w];
    }

    for (const auto &el : BNode.get_forbidden())
        mod_weights.at(el) = std::numeric_limits<dist_type>::max();
    for (const auto &el : BNode.get_required())
        mod_weights.at(el) = -1; //std::numeric_limits<dist_type>::min();


    typedef std::pair<dist_type, int > iPair;
    // Create a priority queue to store vertices that
    // are being preinMST. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // http://geeksquiz.com/implement-min-heap-using-stl/
    std::priority_queue< iPair, std::vector <iPair> , std::greater<iPair> > pq;

    int src = 1; // Taking vertex 0 as source
    size_type n = tsp.size();
    // Create a vector for keys and initialize all
    // keys as infinite (INF)
    std::vector<double> key(n, std::numeric_limits<double>::max()/2.);

    // To store parent array which in turn store MST
    std::vector<int> parent(n, -1);

    // To keep track of vertices included in MST
    std::vector<bool> inMST(n, false);

    // Insert source itself in priority queue and initialize
    // its key as 0.
    pq.push(std::make_pair(0, src));
    key[src] = 0;

    /* Looping till priority queue becomes empty */
    while (!pq.empty()) {
        // The first vertex in pair is the minimum key
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted key (key must be first item
        // in pair)
        int u = pq.top().second;
        pq.pop();

        inMST[u] = true;  // Include vertex in MST

        for (NodeId i = 1 ; i < n ; i++)
        {
            if (i!= u) {
                // Get vertex label and weight of current adjacent
                // of u.
                dist_type weight = mod_weights.at(to_EdgeId(u,i,n));

                //  If v is not in MST and weight of (u,v) is smaller
                // than current key of v
                if (inMST[i] == false && key[i] > weight) {
                    // Updating key of v
                    key[i] = weight;
                    pq.push(std::make_pair(key[i], static_cast<int>(i)));
                    parent[i] = u;
                }
            }
        }
    }

    for (NodeId k = 2; k < n; k++) {
        tree.add_edge(k, parent[k]);
    }
    NodeId smallest = 1;
    for(NodeId k = 2; k < n; k++){
        if(mod_weights.at(to_EdgeId(0, k, n)) < mod_weights.at(to_EdgeId(0, smallest, n))){
            smallest = k;
        }

    }

    NodeId smallest1 =1;
    if(smallest == 1) smallest1 = 2;
    for(NodeId k = 2; k < n; k++){
        if(k != smallest){
            if(mod_weights.at(to_EdgeId(0, k, n)) < mod_weights.at(to_EdgeId(0, smallest1, n))){
                smallest1 = k;
            }
        }
    }
    tree.add_edge(0,smallest);
    tree.add_edge(0,smallest1);
}

template<class coord_type, class dist_type>
dist_type Held_Karp(const TSP::Instance<coord_type, dist_type> &tsp,
                    std::vector<double> & lambda,
                    OneTree &tree,
                    const BranchingNode<coord_type, dist_type> &bn,
                    bool root = false) {
    size_type n = tsp.size();
    std::vector<dist_type> sol_vector;
    std::vector<double> lambda_max(lambda.size(),0), lambda_tmp(lambda) ;

    OneTree tree_max(tree), tree_tmp(tree);
    double t_0 = 0., del_0 = 0., deldel = 0.;
    size_t max_el = 0;
    size_t N = std::ceil(n / 4.) + 5;
    if (root) {
        N =  std::ceil(n*n / 50.) + n +  15;
    }
    compute_minimal_1_tree<coord_type, dist_type>(tree, lambda_tmp, tsp, bn);
    if(root){
        dist_type sum = 0;
        for (const auto &el : tree.get_edges())
            sum += tsp.weight(el);
        t_0 = sum/(2.*n);
    }
    else{
        t_0 = 0;
        for(NodeId i = 0; i < n; i++){
            t_0 += fabs(lambda.at(i));
        }
        t_0 *= 1./(2.*n);
    }
        //t_0 = 1./(2. * n) * std::accumulate(lambda.begin(), lambda.end(), 0.);

    deldel = t_0 / (N*N - N);
    del_0 = 3. * t_0 / (2. * N);

    for (size_t i = 0; i < N; i++) {
        dist_type sum = 0, sum2 = 0;
        for (const auto &el : tree.get_edges())
            sum += tsp.weight(el);
        for (size_t node = 0; node < n; node++)
            sum2 += (tree.get_node(node).degree() - 2.) * lambda_tmp[node];
        sol_vector.push_back(sum + sum2);

        if (i == 0) {
            tree_max = tree;
            lambda_max = lambda_tmp;

            for (size_t j = 0; j < lambda_tmp.size(); j++) {
                lambda_tmp[j] += t_0 * (tree.get_node(j).degree() - 2.);
            }
            t_0 = t_0 - del_0;
            del_0 = del_0 - deldel;
            tree_tmp  = tree;
        }

        if (i > 0) {
            if (sol_vector[max_el] < sol_vector[i]) {
                lambda_max = lambda_tmp;
                max_el = i;
                tree_max = tree;
            }
            for (size_t j = 0; j < lambda_tmp.size(); j++) {
                lambda_tmp[j] += t_0 * (0.6 * (tree.get_node(j).degree() - 2.) + 0.4 * (tree_tmp.get_node(j).degree() - 2.));
            }
            t_0 = t_0 - del_0;
            del_0 = del_0 - deldel;
            tree_tmp = tree;
        }
        tree = OneTree(0, n);
        compute_minimal_1_tree<coord_type, dist_type>(tree, lambda_tmp, tsp, bn);
    }
    if(root){
        lambda = lambda_max;
    }
    tree = tree_max;
    return std::ceil( (1.-EPS) * (*std::max_element(sol_vector.begin(), sol_vector.end())));
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
  ) : required(), forbidden(), lambda(size, 0), tree(0, size), size(tsp.size()), required_neighbors(size),
  forbidden_neighbors(size) {
      HK = Held_Karp(tsp, this->lambda, this->tree, *this, true);
  }

  BranchingNode(const BranchingNode& obj) : tree(0,obj.size) ,       required_neighbors(obj.required_neighbors),
                                            forbidden_neighbors(obj.forbidden_neighbors)
  {
      required = obj.required;
      required_neighbors = obj.required_neighbors;
      forbidden = obj.forbidden;
      size = obj.size;
//      required_neighbors.at(i).neighbors() = obj.required_neighbors.at(i).neighbors();
//      forbidden_neighbors.at(i).neighbors() = obj.forbidden_neighbors.at(i).neighbors();

      tree = obj.tree;
      lambda = obj.lambda;
      HK = obj.HK;
  }

  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1
  ) : required(BNode.get_required()),
      forbidden(BNode.get_forbidden()),
      lambda(BNode.get_lambda()),
      tree(0, tsp.size()),
      size(tsp.size()),
      required_neighbors(BNode.required_neighbors),
      forbidden_neighbors(BNode.forbidden_neighbors)
      {

      add_forbidden(e1);
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }

  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1,
                EdgeId e2
  ) : size(tsp.size()),
      required(BNode.get_required()),
      forbidden(BNode.get_forbidden()),
      lambda(BNode.get_lambda()),
      tree(0, tsp.size()),
      required_neighbors(BNode.required_neighbors),
      forbidden_neighbors(BNode.forbidden_neighbors)
       {
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
      forbidden(BNode.get_forbidden()),
      lambda(BNode.get_lambda()),
      tree(0, tsp.size()),
      required_neighbors(BNode.required_neighbors),
      forbidden_neighbors(BNode.forbidden_neighbors){

      if (both_req) {
          add_required(e1);
          add_required(e2);
      }
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }

  bool operator<( const BranchingNode<coord_type, dist_type>& rhs) const;

  const EdgeId reverse_edge(EdgeId e, size_type n) const {
      NodeId i = 0, j = 0;
      to_NodeId(e, i , j , n);
      return j* n + i;
  }

  bool is_required(EdgeId id) const {
//      assert((std::find(required.begin(), required.end(), id) == required.end() )
//          == ( std::find(required.begin(), required.end(), reverse_edge(id,size)) == required.end()) );
      return (std::find(required.begin(), required.end(), id) != required.end()
      && std::find(required.begin(), required.end(), reverse_edge(id,size)) != required.end())
          ;
  }

  bool is_forbidden(EdgeId id) const {
      return (std::find(forbidden.begin(), forbidden.end(), id) != forbidden.end()
          && std::find(forbidden.begin(), forbidden.end(), reverse_edge(id,size)) != forbidden.end());
  }



  void forbid(NodeId idx, EdgeId e1, EdgeId e2) {
//      EdgeId edge_start = to_EdgeId(idx, 0, size);
      for (NodeId k = 0; k < size; k++) {
          if (idx != k) {
              EdgeId edge = to_EdgeId(idx, k, size);
              if (edge != e1 && edge != e2) {
                  add_forbidden(edge);
              }
          }
      }
  }

  void admit(NodeId idx) {
      for (NodeId k = 0; k < size; k++){
          if (idx != k) {
              EdgeId edge = to_EdgeId(idx, k, size);
              if (!is_forbidden(edge))
                  add_required(edge);
          }
      }
  }


  bool push_required(EdgeId e) {
      if (is_required(e))
          return false;
      assert(!is_forbidden(e));
      NodeId i = 0, j = 0;
      to_NodeId(e, i, j, size);
      required.push_back(e);
      required_neighbors.at(i).add_neighbor(j);
      required.push_back(reverse_edge(e, size));
      required_neighbors.at(j).add_neighbor(i);
      return true;
  }


  void add_required(EdgeId e) {
      // hier pushen wir die edges entsprechend
      if (!push_required(e))
          return;
      NodeId i = 0, j= 0;
      to_NodeId(e,i,j,size);

      if (required_neighbors.at(i).degree() == 2)
          forbid(i,to_EdgeId(i,required_neighbors.at(i).neighbors().at(0),size),
                 to_EdgeId(i,required_neighbors.at(i).neighbors().at(1),size));
      if (required_neighbors.at(j).degree() == 2)
          forbid(j,to_EdgeId(j,required_neighbors.at(j).neighbors().at(0),size),
                 to_EdgeId(j,required_neighbors.at(j).neighbors().at(1),size));
  }

  bool push_forbidden(EdgeId e) {
      if (is_forbidden(e))
          return false;
      assert(!is_required(e));
      NodeId i = 0, j= 0;
      to_NodeId(e,i,j,size);
      forbidden.push_back(e);
      forbidden_neighbors[i].add_neighbor(j);
      forbidden.push_back(reverse_edge(e,size));
      forbidden_neighbors[j].add_neighbor(i);
      return true;
  }


  void add_forbidden(EdgeId e) {
      // hier pushen wir die edges entsprechend
      if( !push_forbidden(e))
          return;
      NodeId i = 0, j= 0;
      to_NodeId(e,i,j,size);

      if (forbidden_neighbors.at(i).degree() == size- 3) {
          admit(i);
      }
      if (forbidden_neighbors.at(j).degree() == size- 3) {
          admit(j);
      }

  }


/*
//  void push_required(EdgeId e, const BranchingNode<coord_type, dist_type> &BNode) {
//      if (is_required(e))
//          return;
//      if (!is_forbidden(e)) {
//          this->push_required(e);
//
//          int count = 0;
//          NodeId i = 0, j = 0;
//          EdgeId e1 = 0;
//          to_NodeId(e,i,j,size);
//          for (NodeId k = 0; k < size; k++){
//              if(k != i && k != j) {
//                  if (is_required(to_EdgeId(i, k, size))) {
//                      count++;
//                      if(count == 1){
//                          e1 = to_EdgeId(i,k,size);
//                      }
//                  }
//              }
//          }
//          if(count == 1){
//              forbid(i, e, e1);
//          }
//          count = 0;
//          for (NodeId k = 0; k < size; k++){
//              if(k!=j && k!= i) {
//                  if (is_required(to_EdgeId(j, k, size))) {
//                      count++;
//                      if(count == 1){
//                          e1 = to_EdgeId(j,k,size);
//                      }
//                  }
//              }
//          }
//          if(count == 1){
//              forbid(j, e, e1);
//              return;
//          }
//
////          to_NodeId(e, i, j, size);
////          for (const auto & el : required) {
////              NodeId idx1 = 0, idx2 = 0;
////              to_NodeId(el,idx1,idx2,size);
////              if (i == idx1) {
////                  forbid(i, e, el);
////              }
////              if (i == idx2) {
////                  forbid(i, e, el);
////              }
////              if (j == idx1) {
////                  forbid(j, e, el);
////              }
////              if (j == idx2) {
////                  forbid(j, e, el);
////              }
////          }
//      }
//  }
//
//  void push_required(EdgeId e) {
//      if (is_required(e))
//          return;
//      if (!is_forbidden(e)) {
//          this->required.push_back(e);
//          this->required.push_back(reverse_edge(e,size));
//      }
//  }
//
//
//    void push_forbidden(EdgeId e, const BranchingNode<coord_type, dist_type> &BNode) {
//        if (is_forbidden(e))
//            return;
//        if (!is_required(e)) {
//            this->push_forbidden(e);
//
//            int count = 0;
//            NodeId i = 0, j = 0;
//            EdgeId e1 = 0, e2 = 0;
//            to_NodeId(e,i,j,size);
//            for (NodeId k = 0; k < size; k++){
//                if(k != i) {
//                    if (!is_forbidden(to_EdgeId(i, k, size))) {
//                        count++;
//                        if(count == 1){
//                            e1 = to_EdgeId(i,k,size);
//                        }
//                        if(count == 2){
//                            e2 = to_EdgeId(i,k,size);
//                        }
//                    }
//                }
//            }
//            if(count == 2){
//                push_required(e1);
//                push_required(e2);
//            }
//            count = 0;
//            for (NodeId k = 0; k < size; k++){
//                if(k!=j) {
//                    if (!is_forbidden(to_EdgeId(j, k, size))) {
//                        count++;
//                        if(count == 1){
//                            e1 = to_EdgeId(j,k,size);
//                        }
//                        if(count == 2){
//                            e2 = to_EdgeId(j,k,size);
//                        }
//                    }
//                }
//            }
//            if(count == 2){
//                push_required(e1);
//                push_required(e2);
//                return;
//            }
//        }
//    }
//
//
//    void push_forbidden(EdgeId e) {
//        if (is_forbidden(e))
//            return;
//        if (!is_required(e)) {
//            forbidden.push_back(e);
//            forbidden.push_back(reverse_edge(e, size));
//        }
//    }
*/
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

  const std::vector<Node> & get_required_neighbors() const {
      return required_neighbors;
  }

  bool tworegular() {
      for (const auto &el : tree.get_nodes())
          if (el.degree() != 2)
              return false;
      return true;
  }

 private:
//  friend class TSP::Instance<coord_type, dist_type>;
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