#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <string>
#include <climits>

namespace TSP {
    //TODO: Adjust constructors

    using size_type = std::size_t;
    using NodeId = size_type ;


    template <class cood_type, class dist_type>
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
                        strstr >> dimnesion;
                    }

                    if (option == "NODE_COORD_SECTION") {
                        scan = true;
                    }
                } while (!scan && file.good());

                if (!scan)
                    throw std::runtime_error("File not in right format");
                std::vector<double> x(dimension), y(dimension);
                while (file.good()) {
                    getline(file, line);
                    std::stringstream strstr;
                    strstr << line;
                    strstr >> option;
                    if (option != "EOF") {
                        size_type id = std::numeric_limits<size_type>;
                        double x = std::numeric_limits<double> , y = std::numeric_limits<double> ;
                        try {
                            strstr >> id >> coord_x >> coord_y;
                            x.push_back(coord_x) , y.push_back(coord_y);
                        }
                        catch (int e) {
                            cout << "An exception occurred while reading the file. Exception Nr. " << e << '\n';
                        }
                    }
                    else break;
                }
                file.close();
                for (size_t i = 0; i < dimension; i++)
                    for (size_t j = 0; j < dimension; j++)
                        if (i != j)
                            _weights.push_back(
                                    distance(x[i],y[i],x[j],y[j])
                            );
            }

            dist_type distance(x1,y1,x2,y2) {
                return std::lround(std::sqrt((_nodes[node1].x - _nodes[node2].x) * (_nodes[node1].x - _nodes[node2].x)
                                             + (nodes[node1].y - _nodes[node2].y) * (_nodes[node1].y - _nodes[node2].y)
                                            )
                                  );
            }

        private:
            std::vector<NodeId> _nodes;
            std::vector<dist_type> _weights;
            size_type dimension;

        };
}