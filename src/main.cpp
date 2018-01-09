#include <iostream>
#include <cstring>
#include <ctime>
#include "../header/tsp.hpp"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "No parameters were given. Please give an --instance ./instance.tsp as an program argument" << std::endl;
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "--instance") != 0) {
        std::cout << "First argument should be an instance of TSP. Execute like ./program --instance ./dir_to_instance.tsp";
        return EXIT_FAILURE;
    }

    std::string file = argv[2];
    std::clock_t begin = clock();
    TSP::Instance<double,double> myTSP(file);
    myTSP.compute_optimal_tour();
    std::clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cerr << "reading, initializing and doing nothing  s" << elapsed_secs << std::endl;
    return EXIT_SUCCESS;
}

//Hello