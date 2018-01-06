#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <climits>
#include <vector>
#include "../src/util.cpp"
#include <sstream>

namespace TSP {
    //TODO: Adjust constructors

    using size_type = std::size_t;
    using NodeId = size_type ;


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
                    if (option != "EOF") {
                        try {
                            strstr >> coord_x >> coord_y;
                            std::cout << option << ' ' << coord_x << ' ' << coord_y << std::endl;
                            _nodes.push_back(std::stoi(option)) , x.push_back(coord_x) , y.push_back(coord_y);
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

        private:
            std::vector<NodeId> _nodes;
            std::vector<dist_type> _weights;
            size_type dimension;

        };
}