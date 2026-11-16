#pragma once

#include "model.hpp"

#include <ostream>

namespace census {

template<class Values>
void write_list(std::ostream& output, const Values& values, std::size_t first, std::size_t end) {
    if (first == end) output << "[]";
    for (std::size_t i = first; i < end; ++i) {
        if (i != first) output << ' ';
        output << values[i];
    }
}

template<class Values>
void write_list(std::ostream& output, const Values& values) {
    write_list(output, values, 0, values.size());
}

void write_fibres(std::ostream& output, const std::vector<LocalFibre>& fibres,
    bool component_cycles = false);

}
