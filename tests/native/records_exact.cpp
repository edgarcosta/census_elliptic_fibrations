#include "records.hpp"

#include <iostream>

int main() {
    using Z = boost::multiprecision::cpp_int;
    using Q = boost::multiprecision::cpp_rational;
    census::InvariantRow row{};
    row.q = 3;
    row.key = UINT64_MAX;
    row.kind = "smooth";
    row.minimality = "everywhere minimal";
    row.status = "partial";
    row.degree = 4;
    row.sign = -1;
    row.lower_l = std::vector<Z>{-((Z(1)<<80)+1), 0};
    row.fibres = std::vector<census::LocalFibre>{
        {"81", 4, "I", 6, 6, 1, 2, -1, {}},
        {"inf", 1, "II*", 0, 12, 4, 1, 0, {}}};
    row.torsion = std::vector<int>{2, 4};
    row.sha_reg = Q(6*((Z(1)<<90)+1))/8;
    row.admissible_sha = std::vector<Z>{1, 4, Z("12157665459056928801")};
    row.split_method = "unresolved";
    row.sha2_dimensions = std::vector<int>{};
    row.message = "trace, \"pending\"\nretry";
    census::CompressedRecords output(std::cout);
    output.write(census::invariant_csv_header()+"\n");
    output.write(census::invariant_csv(row)+"\n");
    row.lower_l = std::vector<Z>{};
    row.degree = 0;
    row.sign = 1;
    row.fibres = std::vector<census::LocalFibre>{};
    row.torsion = std::vector<int>{};
    row.sha_reg = 3;
    row.admissible_sha = std::vector<Z>{};
    row.status = "error";
    row.refinement_status = "error";
    output.write(census::invariant_csv(row)+"\n");
    output.finish();
}
