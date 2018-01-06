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
#include "util.hpp"




namespace TSP {
    using size_type = std::size_t;
    using NodeId = size_type;
    using EdgeId = size_type;

    EdgeId to_EdgeId(NodeId i , NodeId j, size_type N) {
        if ( i == j)
            throw std::runtime_error("Loops are not contained in this instace");
        if (i>j)
            std::swap(i,j);
        return (i)*(N-1) +j -1;
    }

    void to_NodeId(EdgeId e , NodeId & i, NodeId & j, size_type N) {
        i = e/(N-1);
        j = e % (N-1) +1;
    }

    template <class dist_type>
    class BranchingNode;

    template <class coord_type, class dist_type>
    class Instance;

    template <class coord_type, class dist_type>
    void compute_minimal_1_tree(std::vector<EdgeId> & tree,
                                const std::vector<dist_type> lambda,
                                const Instance<double, dist_type> & tsp,
                                const BranchingNode<dist_type> & bn){
        Union_Find uf(tsp.size());
        NodeId i = 0,j = 0 ;
        std::vector<EdgeId> sorted(tsp.num_edges());
        for (size_type index = 0 ; index < tsp.num_edges(); index++)
            sorted[index] = index;
        for (auto const & el : bn.get_required()) {
            tree.push_back(el);
            to_NodeId(el,i,j,tsp.size());
            uf._union(i,j);
        }
        std::stable_sort(sorted.begin(), sorted.end(),
                         [&](EdgeId ind1, EdgeId ind2) {return (tsp.weight(ind1) < tsp.weight(ind2));} );


    return;

    }

    template <class dist_type>
    dist_type Held_Karp(const TSP::Instance<double,dist_type> & tsp,
                        const std::vector<dist_type> & lambda,
                        std::vector<EdgeId> & tree,
                        const BranchingNode<dist_type> & bn) {
        size_type n = tsp.size();
        for (size_t i = 0; i < std::ceil(n*n/50 ) + n +15; i++) {
            compute_minimal_1_tree<double, dist_type>(tree,lambda, tsp, bn);
        }
        return 0;
    }



    template <class coord_type, class dist_type>
    class Instance {
    public:
        Instance(const std::string &filename) {
            std::ifstream file(filename);

            if (!file.is_open())
                throw std::runtime_error("File " + filename + " could not be opened");

            std::string line = "", option = "";
            bool scan = false;
            do {
                getline(file, line);
                stripColons(line);
                std::stringstream strstr;
                strstr << line;
                strstr >> option;

                if (option == "DIMENSION") {
                    strstr >> this->dimension;
                }

                if (option == "NODE_COORD_SECTION") {
                    scan = true;
                }
            } while (!scan && file.good());

            if (!scan)
                throw std::runtime_error("File not in right format");

            std::vector<coord_type> x, y;
            x.reserve(dimension), y.reserve(dimension);
            //_nodes.reserve(dimension);
            //_weights.reserve(dimension*(dimension-1));
            //size_type id = std::numeric_limits<size_type>::max();
            coord_type coord_x = std::numeric_limits<coord_type>::max() , coord_y = std::numeric_limits<coord_type>::max() ;
            while (file.good()) {
                getline(file, line);
                std::stringstream strstr;
                strstr << line;
                strstr >> option;
                // TODO
                if (option != "EOF") {
                    try {
                        strstr >> coord_x >> coord_y;
                        std::cout << option << ' ' << coord_x << ' ' << coord_y << std::endl;
                        _nodes.push_back(std::stoi(option)-1) , x.push_back(coord_x) , y.push_back(coord_y);
                    }
                    catch (int e) {
                        std::cerr << "An exception occurred while reading the file. Exception Nr. " << e << '\n';
                    }
                }
                else break;
            }
            file.close();
            for (size_t i = 0; i < dimension; i++)
                for (size_t j = 0; j < dimension; j++)
                    if (i != j)
                        this->_weights.push_back(
                                distance(x[i], y[i], x[j], y[j])
                        );

        }

        dist_type distance(coord_type x1,coord_type y1,coord_type x2,coord_type y2) {
            return std::lround(std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
        }


        void compute_optimal_tour()  {
            dist_type upperBound = std::numeric_limits<dist_type>::max(); // we want to do a naive tsp tour!
            auto cmp = [](BranchingNode<dist_type>, BranchingNode<dist_type>) { return true; };
            std::priority_queue<BranchingNode<dist_type>,
                    std::vector<BranchingNode<dist_type> >, decltype(cmp)> Q(cmp);
            Q.push(BranchingNode<dist_type>(std::vector<EdgeId>(), std::vector<EdgeId>(),*this,std::vector<dist_type>(dimension)));

            while (!Q.empty()) {
                BranchingNode<dist_type> current_BN(Q.top());
                Q.pop();
                if (current_BN.HK > upperBound)
                    continue;
                else {

                }
            }
        }

        size_type size() const {
            return  dimension;
        }

        size_type num_edges() const {
            return _weights.size();
        }

        dist_type weight(EdgeId id) const {
            return _weights[id];
        }

    private:
        std::vector<NodeId> _nodes;
        std::vector<dist_type> _weights;
        size_type dimension;

    };

    template <class dist_type>
    class BranchingNode {
    public:
        BranchingNode(const std::vector<EdgeId > & req,
                      const std::vector<EdgeId> & forb,
                      const Instance<double, dist_type> & tsp,
                      const std::vector<dist_type> & lambda
        ) : required(req) , forbidden(forb), lambda(lambda) {

            this->HK = Held_Karp(tsp, lambda, tree, *this);
        }

        bool check_required(size_type id);
        bool check_forbidden(size_type id);

      const std::vector<EdgeId> & get_required() const {
          return required;
      }
      const std::vector<EdgeId> & get_forbidden() const {
          return forbidden;
      }




    private:
        friend Instance<double,dist_type>;
        std::vector<EdgeId> required;
        std::vector<EdgeId> forbidden;
        std::vector<dist_type> lambda;
        std::vector<EdgeId> tree;
        dist_type HK;
    };
}


#endif // BRANCHANDBOUNDTSP_TSP_HPP