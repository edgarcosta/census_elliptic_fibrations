#pragma once

#include "../invariants/invariants.hpp"

#include <array>
#include <optional>
#include <ostream>

namespace census {

std::string invariant_csv_header();
std::string invariant_csv(const Model& model, const InvariantData& data,
    const std::optional<std::string>& error);

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
