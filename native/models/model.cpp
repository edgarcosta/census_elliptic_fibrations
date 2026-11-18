#include "model.hpp"
#include "weierstrass.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace census {
namespace {
using NTL::zz_pX;
bool form_power(const zz_pX& f, int exponent) {
    for (long i = 0; i <= NTL::deg(f); ++i) {
        if (i % exponent != 0 && !NTL::IsZero(NTL::coeff(f, i))) return false;
    }
    return true;
}

std::string fibre_kind(const Coefficients& a, unsigned q) {
    const auto& [a1, a2, a3, a4, a6] = a;
    const zz_pX discriminant = WeierstrassInvariants(a).discriminant();
    bool cusp;
    bool regular;
    if (q == 2) {
        cusp = NTL::IsZero(a1) && NTL::IsZero(a3);
        regular = !form_power(a4, 2) || !form_power(a6, 2);
    } else {
        cusp = NTL::IsZero(a2) && NTL::IsZero(a4);
        regular = !form_power(a6, 3);
    }
    if (!NTL::IsZero(discriminant)) return "smooth";
    return cusp && regular ? "regular cusp" : "nonregular";
}

std::vector<zz_pX> residues(unsigned q, int width) {
    unsigned count = 1;
    for (int i = 0; i < width; ++i) count *= q;
    std::vector<zz_pX> result;
    result.reserve(count);
    for (unsigned code = 0; code < count; ++code) {
        zz_pX f;
        unsigned digits = code;
        for (int i = 0; i < width; ++i) {
            NTL::SetCoeff(f, i, digits % q);
            digits /= q;
        }
        result.push_back(std::move(f));
    }
    return result;
}

bool divisible(const zz_pX& f, const zz_pX& modulus) {
    return NTL::IsZero(f % modulus);
}

bool nonminimal_at(const Coefficients& a, const zz_pX& pi, unsigned q) {
    const auto& [a1, a2, a3, a4, a6] = a;
    const int degree = NTL::deg(pi);
    const zz_pX pi2 = NTL::sqr(pi);
    const zz_pX pi4 = NTL::sqr(pi2);
    const zz_pX pi6 = pi4*pi2;
    if (q == 2) {
        if (!divisible(a1, pi)) return false;
        const zz_pX pi3 = pi2*pi;
        const auto shifts = residues(q, 3*degree);
        for (const auto& slope : residues(q, degree)) {
            const zz_pX r = slope*a1 + NTL::sqr(slope);
            if (!divisible(a3 + r*a1, pi3)) continue;
            for (const auto& v : shifts) {
                const zz_pX c4 = slope*a3 + (v + r*slope)*a1 + NTL::sqr(r);
                if (!divisible(a4 + c4, pi4)) continue;
                const zz_pX c6 = NTL::power(r, 3) + v*a3 + NTL::sqr(v) + r*v*a1;
                if (divisible(a6 + r*a4 + c6, pi6)) return true;
            }
        }
    } else {
        if (!divisible(a2, pi2)) return false;
        for (const auto& r : residues(q, 2*degree)) {
            if (!divisible(a4 - a2*r, pi4)) continue;
            if (divisible(a6 + a4*r + a2*NTL::sqr(r) + NTL::power(r, 3), pi6)) {
                return true;
            }
        }
    }
    return false;
}

bool everywhere_minimal(const Coefficients& affine, const Coefficients& infinity,
                        unsigned q) {
    const bool cusp = q == 2
        ? NTL::IsZero(affine[0]) && NTL::IsZero(affine[2])
        : NTL::IsZero(affine[1]) && NTL::IsZero(affine[3]);
    if (cusp && (q == 2
        ? form_power(affine[3], 4) && form_power(affine[4], 2)
        : form_power(affine[4], 3))) return false;

    const zz_pX variable = polynomial({0, 1});
    if (nonminimal_at(infinity, variable, q)) return false;
    for (unsigned c = 0; c < q; ++c) {
        if (nonminimal_at(affine, polynomial({static_cast<int>(c), 1}), q)) return false;
    }
    // The degree bounds in the homogeneous local conditions exclude higher places.
    if (!cusp) {
        const std::vector<std::vector<int>> quadratics = q == 2
            ? std::vector<std::vector<int>>{{1, 1, 1}}
            : std::vector<std::vector<int>>{{1, 0, 1}, {2, 1, 1}, {2, 2, 1}};
        for (const auto& coefficients : quadratics) {
            if (nonminimal_at(affine, polynomial(coefficients), q)) return false;
        }
    }
    return true;
}

void trim(std::vector<int>& coefficients) {
    while (coefficients.size() > 1 && coefficients.back() == 0) coefficients.pop_back();
}
}

Model model_from_key(unsigned q, std::uint64_t key) {
    if (q != 2 && q != 3) throw std::invalid_argument("q must be the prime 2 or 3");
    const std::uint64_t bound = q == 2 ? 4294967296ULL : 7625597484987ULL;
    if (key >= bound) throw std::invalid_argument("key is outside the encoding range");

    NTL::zz_pPush context(q);
    Model model;
    model.q = q;
    model.key = key;
    const std::array<bool, 5> included = q == 2
        ? std::array<bool, 5>{true, false, true, true, true}
        : std::array<bool, 5>{false, true, false, true, true};
    for (int i = 4; i >= 0; --i) {
        model.affine[i].resize(model_degrees[i] + 1);
        if (!included[i]) continue;
        for (int& digit : model.affine[i]) {
            digit = key % q;
            key /= q;
        }
    }
    for (int i = 0; i < 5; ++i) {
        model.infinity[i] = model.affine[i];
        std::reverse(model.infinity[i].begin(), model.infinity[i].end());
        trim(model.affine[i]);
        trim(model.infinity[i]);
    }
    const auto affine = polynomials(model.affine);
    const auto infinity = polynomials(model.infinity);
    model.kind = fibre_kind(affine, q);
    model.minimality = everywhere_minimal(affine, infinity, q)
        ? "everywhere minimal" : "nonminimal somewhere";
    model.status = model.kind == "smooth" && model.minimality == "everywhere minimal"
        ? "partial" : "unsupported";
    return model;
}

}
