#pragma once

#include <NTL/lzz_pX.h>

#include <array>
#include <vector>

namespace census {

using Coefficients = std::array<NTL::zz_pX, 5>;

struct WeierstrassInvariants {
    NTL::zz_pX b2, b4, b6, b8;

    explicit WeierstrassInvariants(const Coefficients& a) {
        const auto& [a1, a2, a3, a4, a6] = a;
        b2 = a1*a1 + 4*a2;
        b4 = a1*a3 + 2*a4;
        b6 = a3*a3 + 4*a6;
        b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4;
    }

    NTL::zz_pX discriminant() const {
        return -b2*b2*b8 - 8*b4*b4*b4 - 27*b6*b6 + 9*b2*b4*b6;
    }
};

inline NTL::zz_pX polynomial(const std::vector<int>& coefficients) {
    NTL::zz_pX result;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        NTL::SetCoeff(result, i, coefficients[i]);
    }
    return result;
}

inline Coefficients polynomials(const std::array<std::vector<int>, 5>& coefficients) {
    Coefficients result;
    for (int i = 0; i < 5; ++i) result[i] = polynomial(coefficients[i]);
    return result;
}

}
