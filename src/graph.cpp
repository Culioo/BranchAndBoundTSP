#include "../header/graph.hpp"

#include <ostream>
#include <stdexcept>

namespace TSP {
/////////////////////////////////////////////
//! \c Node definitions
/////////////////////////////////////////////

    void Node::add_neighbor(NodeId const id) {
        _neighbors.push_back(id);
    }


} // namespace TSP
