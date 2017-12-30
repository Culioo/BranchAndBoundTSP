#include <iostream>
#include <cstring>
#include <ctime>


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "No parameters were given. Please give an --instance ./instance.tsp as an program argument" << std::endl;
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "--instance") != 0) {
        std::cout << "First argument should be an instance of TSP. Execute like ./program --instance ./dir_to_instance.tsp";
        return EXIT_FAILURE;
    }


    std::clock_t begin = clock();

    std::clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "reading, initializing and computing the matching took  s" << elapsed_secs << std::endl;
    return EXIT_SUCCESS;
}

