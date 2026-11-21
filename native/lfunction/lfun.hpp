#pragma once

#include "field.hpp"

#include <vector>

namespace census {

struct TraceTable {
    TraceField field;
    std::vector<int> ordinary;
    int nodal;

    TraceTable(int p, int degree);
};

}
