#pragma once

#include <cstdint>
#include <vector>

namespace census {

struct TraceField {
    int p, degree, size, order, mask;
    std::vector<int> polynomial, logarithms, traces;
    std::vector<std::uint32_t> exponents, index_mask;

    TraceField(int characteristic, int extension_degree);
    std::uint32_t add(std::uint32_t a, std::uint32_t b) const;
    std::uint32_t negate(std::uint32_t a) const;
    int index(std::uint32_t a) const;
    std::uint32_t exp(std::int64_t e) const;
    std::uint32_t multiply(std::uint32_t a, std::uint32_t b) const;
    int trace(std::uint32_t a) const { return traces[index(a)]; }
};

}
