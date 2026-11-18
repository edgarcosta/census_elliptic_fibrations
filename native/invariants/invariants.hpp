#pragma once

#include "../lfunction/lfun.hpp"

#include <iosfwd>

namespace census {

struct TorsionData {
    std::vector<int> factors;
    std::vector<std::string> points;
};

TorsionData complete_torsion(const Model& model, const LocalModels& local);

struct DerivedInvariants {
    int rank = 0, geometric_rank = 0, picard_number, geometric_picard_number;
    std::vector<mpz_class> p2, point_counts;
    int zeta_sign, brauer_height = 0;
    mpz_class artin_tate_product, ns_splitting_degree = 1;
    int j_degree, j_inseparability_exponent = 0;
    mpq_class sha_reg;
};

DerivedInvariants derived_invariants(const Model& model, const LocalModels& local,
    const TorsionData& torsion, const LPolynomial& lfun);

void write_derived_invariants(std::ostream& output, const DerivedInvariants& result);

struct CandidateInvariants {
    std::vector<mpz_class> sha;
    std::vector<mpq_class> regulators;

    std::string_view split_method() const { return sha.size() == 1 ? "constraints" : "unresolved"; }
};

// Supplied derived values and fibres describe the same elliptic K3 surface.
CandidateInvariants admissible_candidates(const DerivedInvariants& data,
    const std::vector<LocalFibre>& fibres);

void write_candidate_invariants(std::ostream& output, const CandidateInvariants& candidates);

}
