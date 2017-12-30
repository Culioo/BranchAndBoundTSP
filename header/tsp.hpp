#include <cstdlib>
#include <fstream>
#include <algorithm>
namespace TSP {
    using size_type = std::size_t;

    using NodeId = size_type ;

    template <class cood_type, class dist_type>
        class Instance {
        public:
            class Node {
            public:
            private:
                friend Instance;
                cood_type x,y;
            };

            Instance() {

            }

            dist_type distance(NodeId node1, NodeId node2) {
                return std::lround(std::sqrt((_nodes[node1].x - _nodes[node2].x) * (_nodes[node1].x - _nodes[node2].x)
                                             + (nodes[node1].y - _nodes[node2].y) * (_nodes[node1].y - _nodes[node2].y)
                                            )
                                  );
            }

        private:
            std::vector<NodeId> _nodes;
        };
    template <class cood_type, class dist_type>
    Instance<cood_type, dist_type> read_instance(const std::string &filename) {
        std::ifstream file(filename);

        if (!file.is_open())
            throw std::runtime_error("File " + filename + " could not be opened");

        NodeId num_nodes = 0;
        std::string line = "", option = "";
        bool done = false;
        do {
            getline(file, line);
            std::stringstream strstr;
            strstr << line;
            strstr >> option;

            if (option == "p") {
                strstr >> option;
                if (option == "edge") {

                    strstr >> num_nodes >> num_edges;

                    done = true;
                }
            }
        } while (!done && file.good());

        if (!done)
            throw std::runtime_error("File not in right format");

        ED::Graph result{static_cast<ED::NodeId >(num_nodes)};
        while (file.good()) {
            getline(file, line);
            std::stringstream strstr;
            strstr << line;
            strstr >> option;
            if (option != "c") {
                if (option == "e") {
                    ED::DimacsId v = ED::invalid_dimacs_id, w = ED::invalid_dimacs_id;
                    strstr >> v >> w;
                    if (v == ED::invalid_dimacs_id || w == ED::invalid_dimacs_id)
                        throw std::runtime_error("One node ID was out of bounds");
                    result.add_edge(ED::from_dimacs_id(v), ED::from_dimacs_id(w));
                }
            }
            option = "";
        }
        file.close();

        return result;
    }
    };


}