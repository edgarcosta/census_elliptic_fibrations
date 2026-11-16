#include "format.hpp"

namespace census {

void write_fibres(std::ostream& output, const std::vector<LocalFibre>& fibres,
    bool component_cycles) {
    if (fibres.empty()) output << "[]";
    for (std::size_t i = 0; i < fibres.size(); ++i) {
        if (i) output << ';';
        const auto& f = fibres[i];
        output << f.place << ':' << f.degree << ':' << f.type << ':' << f.index
               << ':' << f.discriminant_valuation << ':' << f.conductor << ':'
               << f.tamagawa << ':' << f.trace;
        if (component_cycles) {
            output << ':';
            write_list(output, f.component_cycles);
        }
    }
}

}
