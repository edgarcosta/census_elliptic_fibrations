#pragma once

#include "model.hpp"

#include <ostream>

namespace census {

template<class Values>
void write_joined(std::ostream& output, const Values& values, char separator,
    std::size_t first, std::size_t end) {
    for (std::size_t i = first; i < end; ++i) {
        if (i != first) output << separator;
        output << values[i];
    }
}

template<class Values>
void write_joined(std::ostream& output, const Values& values, char separator) {
    write_joined(output, values, separator, 0, values.size());
}

template<class Values>
void write_list(std::ostream& output, const Values& values, std::size_t first, std::size_t end) {
    if (first == end) output << "[]";
    else write_joined(output, values, ' ', first, end);
}

template<class Values>
void write_list(std::ostream& output, const Values& values) {
    write_list(output, values, 0, values.size());
}

void write_fibres(std::ostream& output, const std::vector<LocalFibre>& fibres,
    bool component_cycles = false);

}
