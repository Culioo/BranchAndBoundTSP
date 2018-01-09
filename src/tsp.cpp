#include "../header/tsp.hpp"
#include <fstream>

namespace TSP{
    void Instance<double ,double >::print_optimal_length() {
        int length = 0;
        for (int count = 0; count < this->size() - 1; count++) {
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

    void Instance<double ,double >::print_optimal_tour(const std::string &filename){

        std::ofstream file_to_print;
        file_to_print.open(filename, std::ios::out);


        file_to_print << "TYPE : TOUR" << std::endl;
        file_to_print << "DIMENSION : " << this->size() << std::endl;
        file_to_print << "TOUR_SECTION" << std::endl;
        for( int i = 0; i < this->size(); i++){
            file_to_print << this->_tour.at(i) << std::endl;
        }
        file_to_print << "-1" << std::endl;
        file_to_print << "EOF" << std::endl;
    }
}
