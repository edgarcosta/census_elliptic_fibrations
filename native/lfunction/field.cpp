#include "field.hpp"

#include <stdexcept>

namespace census {
namespace {
constexpr int polynomials2[12][13] = {
    {1,1}, {1,1,1}, {1,1,0,1}, {1,1,0,0,1}, {1,0,1,0,0,1},
    {1,1,0,1,1,0,1}, {1,1,0,0,0,0,0,1}, {1,0,1,1,1,0,0,0,1},
    {1,0,0,0,1,0,0,0,0,1}, {1,1,1,1,0,1,1,0,0,0,1},
    {1,0,1,0,0,0,0,0,0,0,0,1}, {1,1,0,1,0,1,1,1,0,0,0,0,1}
};
constexpr int polynomials3[12][13] = {
    {1,1}, {2,2,1}, {1,2,0,1}, {2,0,0,2,1}, {1,2,0,0,0,1},
    {2,2,1,0,2,0,1}, {1,0,2,0,0,0,0,1}, {2,2,2,0,1,2,0,0,1},
    {1,1,2,2,0,0,0,0,0,1}, {2,1,0,0,2,2,2,0,0,0,1},
    {1,0,2,0,0,0,0,0,0,0,0,1}, {2,0,1,0,1,1,1,0,0,0,0,0,1}
};
}

std::uint32_t TraceField::add(std::uint32_t a, std::uint32_t b) const {
    if (p == 2) return a ^ b;
    const auto al = a & 65535, ah = a >> 16, bl = b & 65535, bh = b >> 16;
    const auto an = al | ah, bn = bl | bh;
    return ((al & ~bn) | (bl & ~an) | (ah & bh))
        | (((ah & ~bn) | (bh & ~an) | (al & bl)) << 16);
}

std::uint32_t TraceField::negate(std::uint32_t a) const {
    return p == 2 ? a : (a >> 16) | (a << 16);
}

int TraceField::index(std::uint32_t a) const {
    return p == 2 ? int(a) : int(index_mask[a & mask] + 2 * index_mask[a >> 16]);
}

std::uint32_t TraceField::exp(std::int64_t e) const {
    e %= order;
    if (e < 0) e += order;
    return exponents[e];
}

std::uint32_t TraceField::multiply(std::uint32_t a, std::uint32_t b) const {
    return a && b ? exp(logarithms[index(a)] + logarithms[index(b)]) : 0;
}

TraceField::TraceField(int characteristic, int extension_degree)
    : p(characteristic), degree(extension_degree), size(1) {
    if ((p != 2 && p != 3) || degree < 1 || degree > 12) {
        throw std::invalid_argument("p must be 2 or 3 and degree must be between 1 and 12");
    }
    mask = (1 << degree) - 1;
    for (int i = 0; i < degree; ++i) size *= p;
    order = size - 1;
    index_mask.resize(1 << degree);
    for (int i = 1; i < (1 << degree); ++i) {
        const int bit = __builtin_ctz(unsigned(i));
        int weight = 1;
        for (int j = 0; j < bit; ++j) weight *= 3;
        index_mask[i] = index_mask[i & (i - 1)] + weight;
    }
    const int* coefficients = p == 2 ? polynomials2[degree - 1] : polynomials3[degree - 1];
    std::uint32_t modulus = 0;
    for (int i = 0; i < degree; ++i) {
        if (coefficients[i]) {
            modulus |= std::uint32_t(1) << (i + (p == 3 && coefficients[i] == 2 ? 16 : 0));
        }
    }
    exponents.resize(order);
    logarithms.assign(size, -1);
    std::uint32_t a = 1;
    for (int e = 0; e < order; ++e) {
        exponents[e] = a;
        if (!a || logarithms[index(a)] != -1) {
            throw std::runtime_error("nonprimitive field polynomial");
        }
        logarithms[index(a)] = e;
        const int top = ((a >> (degree - 1)) & 1)
            + (p == 3 ? 2 * ((a >> (degree + 15)) & 1) : 0);
        a = p == 2 ? (a << 1) & mask
                   : ((a << 1) & mask) | (((a >> 16 << 1) & mask) << 16);
        if (top) a = add(a, p == 2 || top == 2 ? modulus : negate(modulus));
    }
    if (a != 1) throw std::runtime_error("field cycle did not close");
    std::vector<int> basis_traces(degree);
    for (int i = 0; i < degree; ++i) {
        std::uint32_t b = 1u << i, t = 0;
        for (int j = 0; j < degree; ++j) {
            t = add(t, b);
            b = exp(std::int64_t(logarithms[index(b)]) * p);
        }
        if (index(t) >= p) throw std::runtime_error("field trace is not constant");
        basis_traces[i] = index(t);
    }
    traces.assign(size, 0);
    for (const auto b : exponents) {
        int t = 0;
        for (int i = 0; i < degree; ++i) {
            t += basis_traces[i] * (((b >> i) & 1)
                + (p == 3 ? 2 * ((b >> (i + 16)) & 1) : 0));
        }
        traces[index(b)] = t % p;
    }
}

}
