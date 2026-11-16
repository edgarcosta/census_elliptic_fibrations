#pragma once

#include "lfun.hpp"

#include <iosfwd>

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

void write_derived_invariants(std::ostream& output, const DerivedInvariants& result);

struct CandidateInvariants {
    std::vector<boost::multiprecision::cpp_int> sha;
    std::vector<boost::multiprecision::cpp_rational> regulators;
};

// Supplied derived values and fibres describe the same elliptic K3 surface.
CandidateInvariants admissible_candidates(const DerivedInvariants& data,
    const std::vector<LocalFibre>& fibres);

void write_candidate_invariants(std::ostream& output, const DerivedInvariants& data,
    const CandidateInvariants& candidates);

}
