#pragma once

#include "lfun.hpp"

namespace census {

struct TorsionData {
    std::vector<int> factors;
    std::vector<std::string> points;
    std::vector<int> orders;
};

TorsionData complete_torsion(const Model& model, const LocalModels& local);

struct DerivedInvariants {
    int rank = 0, geometric_rank = 0, picard_number, geometric_picard_number;
    std::vector<boost::multiprecision::cpp_int> p2, point_counts;
    int zeta_sign, brauer_height = 0;
    boost::multiprecision::cpp_int artin_tate_product, ns_splitting_degree = 1;
    int j_degree, j_inseparability_exponent = 0;
    boost::multiprecision::cpp_rational sha_reg;
};

// Supplied arithmetic is complete over GF(model.q^constant_degree), as in Magma.
DerivedInvariants derived_invariants(const Model& model, const LocalModels& local,
    const TorsionData& torsion, const LPolynomial& lfun, int constant_degree = 1);

}
