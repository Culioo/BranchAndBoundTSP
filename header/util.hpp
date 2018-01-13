//
// Created by leon on 06.01.18.
//

#ifndef BRANCHANDBOUNDTSP_UTIL_HPP
#define BRANCHANDBOUNDTSP_UTIL_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

class Union_Find {
    friend std::ostream& operator<<(std::ostream &os, const Union_Find& n);
public:
    /** Constructor of the Union Find Object with \f$ n \f$ objects
     * where all objects are first initializied with themselves.
     *
     * @param number_of_elements number of elements in the union find structure
     */
    Union_Find(size_t number_of_elements) :
            parent(std::vector<size_t>(number_of_elements)),
            rank(std::vector<size_t>(number_of_elements)),
            n(number_of_elements)
    {
        for (size_t id = 0; id<number_of_elements; id++) {
            parent[id] = id;
            rank[id] = 0;
        }
    }

    /**
     * Unifies the Partition classes of two given representatives
     * @param rep1 First representative
     * @param rep2 Second representative
     */
    void _union(const size_t x, const size_t y) {
        //std::cout << rep1 << ' ' << rep2 << std::endl;
        size_t xroot = _find(x);
        size_t yroot = _find(y);

        if(rank.at(xroot) < rank.at(yroot)){
            parent.at(xroot) = yroot;
        }
        else if(rank.at(xroot) > rank.at(yroot)){
            parent.at(yroot) = xroot;
        }
        else
        {
            parent.at(yroot) = xroot;
            rank.at(xroot)++;
        }
    }
    /**
     * finds the corresopnding partition class to an element
     * @param el
     * @return Representative of el
     */
    size_t _find(const size_t el) {
        if(parent.at(el) != el){
            parent.at(el) = _find(parent.at(el));
        }
        return parent.at(el);

    }

    /**
     * Simple print function to display the current values of the union find structure.
     * This function was primarily used for debugging purposes.
     */
    void print() {
        std::cout << "Union find assignment" << "\n";
        for (size_t id = 0; id < n; id++) {
            std::cout << id << ' ';
        }
        std::cout << '\n';
        for (size_t id = 0; id < n; id++) {
            std::cout << _find(id) << ' ';
        }
        std::cout << '\n';
    }

private:
    std::vector<size_t> representative; //partitioning class in the union_find tree
    std::vector<size_t> parent; //parent in the union_find tree
    std::vector<size_t> rank;//rank in the union_find tree
    size_t n; //number of elements
    //rank in the union_find tree


};


void stripColons(std::string &x) {
    auto it = std::remove_if(std::begin(x),std::end(x),[](char c){return (c == ':');});
    x.erase(it, std::end(x));
}

#endif //BRANCHANDBOUNDTSP_UTIL_HPP

