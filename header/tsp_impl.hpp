//
// Created by adyck on 1/13/18.
//

#ifndef BRANCHANDBOUNDTSP_TSP_IMPL_HPP
#define BRANCHANDBOUNDTSP_TSP_IMPL_HPP

#include <cassert>
namespace TSP {
template<class coord_type, class dist_type>
Instance<coord_type, dist_type>::Instance(const std::string &filename) {
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
    coord_type coord_x = std::numeric_limits<coord_type>::max(), coord_y = std::numeric_limits<coord_type>::max();
    while (file.good()) {
        getline(file, line);
        std::stringstream strstr;
        strstr << line;
        strstr >> option;
        // TODO: Capture if EOF is missing
        if (option != "EOF") {
            try {
                strstr >> coord_x >> coord_y;
                std::cout << option << ' ' << coord_x << ' ' << coord_y << std::endl;
                _nodes.push_back(std::stoi(option) - 1), x.push_back(coord_x), y.push_back(coord_y);
            }
            catch (int e) {
                std::cerr << "An exception occurred while reading the file. Exception Nr. " << e << '\n';
            }
        } else break;
    }
    file.close();
    for (size_t i = 0; i < dimension; i++)
        for (size_t j = 0; j < dimension; j++)
            this->_weights.push_back(
                distance(x[i], y[i], x[j], y[j])
            );
    _tour = std::vector<NodeId >(dimension);
}


template<class coord_type, class dist_type>
void Instance<coord_type, dist_type>::compute_optimal_tour() {
    typedef BranchingNode <coord_type, dist_type> BNode;

    dist_type upperBound = std::numeric_limits<dist_type>::max(); // we want to do a naive tsp tour!
    auto cmp = [](BNode, BNode) { return true; };
    std::priority_queue<BNode,
                        std::vector<BNode>, decltype(cmp)> Q(cmp);
    Q.push(BranchingNode<coord_type, dist_type>(*this)); // Adding empty node to Q

    while (!Q.empty()) {
        BNode current_BNode(Q.top());
        Q.pop();
        if (current_BNode.get_HK() > upperBound)
            continue;
        else {
            if (current_BNode.tworegular()) {
                upperBound = current_BNode.get_HK();
                std::cerr << "Upper Bound " << upperBound << std::endl;
//                _tour = current_BNode.get_tree();
                continue;
            } else {
//                std::cout << "vindaloop" << std::endl;
                size_type gl_i = 0, choice1 = std::numeric_limits<size_type>::max(),
                    choice2 = std::numeric_limits<size_type>::max();
                for (NodeId node = 1;  node <  current_BNode.get_tree().get_nodes().size(); node++) {
                    if (current_BNode.get_tree().get_node(node).degree() > 2) {
                        gl_i = node;
                        break;
                    }
                }

                assert(gl_i!=0);
                size_t counter = 0;
                for (const auto & el : current_BNode.get_tree().get_node(gl_i).neighbors()){
                    if (!current_BNode.is_required(to_EdgeId(gl_i, el, size())) &&
                        !current_BNode.is_forbidden(to_EdgeId(gl_i, el , size()))) {
                        if (counter == 0)
                            choice1 = el;
                        if (counter == 1)
                            choice2 = el;
                        counter++;
                        if (counter > 1)
                            break;
                    }
                }
                BNode q1(current_BNode, *this, to_EdgeId(gl_i, choice1, this->size())),
                    q2(current_BNode,
                       *this,
                       to_EdgeId(gl_i, choice1, this->size()),
                       to_EdgeId(gl_i, choice2, this->size()));
                Q.push(q1);
                Q.push(q2);

                bool take_q3 = true;
                for (auto node_it =  current_BNode.get_tree().get_node(gl_i).neighbors().begin();
                     node_it != current_BNode.get_tree().get_node(gl_i).neighbors().end(); node_it++) {
                    if (current_BNode.is_required(to_EdgeId(gl_i, *node_it, this->size()))) {
                        take_q3 = false;
                        break;
                    }

                }
                if (take_q3) {
                    BNode q3(current_BNode, *this,
                             to_EdgeId(gl_i, choice1, this->size()),
                             to_EdgeId(gl_i, choice2, this->size()),
                             true);
                    Q.push(q3);
                }

                }
            }
        }
    }


template<class coord_type, class dist_type>
void Instance<coord_type, dist_type>::print_optimal_length() {
    int length = 0;
    for (int count = 0; count < this->_tour.size() - 1; count++) {
        NodeId i = this->_tour.at(count);
        NodeId j = this->_tour.at(count + 1);
        EdgeId edge = i * this->size() + j;
        length += this->_weights.at(edge);
    }
    NodeId i = this->_tour.at(this->size() - 1);
    NodeId j = this->_tour.at(0);
    EdgeId edge = i * this->size() + j;
    length += this->_weights.at(edge);

    std::cout << "The found tour is of length " << length << std::endl;
}

template<class coord_type, class dist_type>
void Instance<coord_type, dist_type>::print_optimal_tour(const std::string &filename) {
    if (this->_tour.size() != this->size())
        throw std::runtime_error("No tour computed yet!");

    std::ofstream file_to_print;
    file_to_print.open(filename, std::ios::out);

    file_to_print << "TYPE : TOUR" << std::endl;
    file_to_print << "DIMENSION : " << this->size() << std::endl;
    file_to_print << "TOUR_SECTION" << std::endl;
    std::vector<NodeId> path(size());
    for (size_t i = 0; i < this->size(); i++) {
        NodeId id1 = std::numeric_limits<NodeId>::max(),
            id2 = std::numeric_limits<NodeId>::max();
        to_NodeId(this->_tour.at(i), id1, id2, size());
        path.push_back(id1);
        path.push_back(id2);
    }
    std::stable_sort(path.begin(),path.end());
    auto last =  std::unique(path.begin(),path.end());
    path.erase(last,path.end());
    for (auto it = path.begin(); it != path.end(); it++) {
        file_to_print << *it << std::endl;
    }
    file_to_print << "-1" << std::endl;
    file_to_print << "EOF" << std::endl;
}


} //end namespace TSP

#endif //BRANCHANDBOUNDTSP_TSP_IMPL_HPP
