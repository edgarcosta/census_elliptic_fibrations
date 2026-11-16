#pragma once

#include "invariants.hpp"

#include <array>
#include <optional>
#include <ostream>

namespace census {

struct InvariantRow {
    unsigned q;
    std::uint64_t key;
    std::string kind, minimality, status;
    std::optional<int> degree, sign;
    std::optional<std::vector<boost::multiprecision::cpp_int>> lower_l;
    std::optional<std::vector<LocalFibre>> fibres;
    std::optional<std::vector<int>> torsion;
    std::optional<boost::multiprecision::cpp_rational> sha_reg;
    std::optional<std::vector<boost::multiprecision::cpp_int>> admissible_sha;
    std::string split_method;
    std::optional<std::vector<int>> sha2_dimensions;
    std::string refinement_status, message;
};

InvariantRow invariant_row(const Model& model, LPolynomialCounter& counter);
std::string invariant_csv_header();
std::string invariant_csv(const InvariantRow& row);

class CompressedRecords {
public:
    explicit CompressedRecords(std::ostream& output);
    ~CompressedRecords();
    void write(std::string_view text);
    void finish();
    CompressedRecords(const CompressedRecords&) = delete;
    CompressedRecords& operator=(const CompressedRecords&) = delete;

private:
    struct State;
    std::unique_ptr<State> state;
};

}
