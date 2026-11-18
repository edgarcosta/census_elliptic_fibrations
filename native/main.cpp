#include "lfunction/lfun.hpp"
#include "models/model.hpp"
#include "invariants/invariants.hpp"
#include "records/records.hpp"
#include "records/format.hpp"

#include <charconv>
#include <cstdint>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
void print_help(std::ostream& output) {
    output << "Usage: census-native <command>\n\n"
              "Commands:\n"
              "  help        Show this help\n"
              "  invariants  Compute exact K3 invariants\n"
              "  count       Count points on the ordinary cubic families\n\n"
              "Options:\n"
              "  -h, --help                Show this help\n"
              "  --backend auto|scalar|avx2  Counting backend (default auto)\n"
              "  --diagnostic              Emit expanded uncompressed comparison CSV\n"
              "Input: one q key pair per line, q=2 or 3, decimal uint64 key.\n"
              "Output: version 1 ASCII CSV compressed with zstd level 9.\n"
              "Smooth everywhere minimal rows receive K3 arithmetic; other cells are unsupported.\n";
    output << "\ncount reads p degree parameter rows.\n"
              "Parameter is an exponent 0 <= e < p^degree-1, or zero for parameter 0.\n"
              "The parameter at exponent e is alpha^e for the retained field polynomial.\n"
              "Degrees 1 through 12 are supported. Output: p,degree,parameter,trace,points.\n"
              "Equations: y^2+x*y=x^3+a over F_2; y^2=x^3+x^2+a over F_3.\n";
}

std::uint64_t unsigned_integer(std::string_view value, const char* field) {
    std::uint64_t result;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size()) {
        throw std::invalid_argument(std::string(field) + " must be a decimal uint64 integer");
    }
    return result;
}

void write_coefficients(const std::array<std::vector<int>, 5>& chart) {
    for (const auto& coefficients : chart) {
        std::cout << ',';
        census::write_joined(std::cout, coefficients, ' ');
    }
}

// Completed stages remain available if a later computation throws.
void compute_invariants(const census::Model& model, census::LPolynomialCounter& counter, census::InvariantData& data) {
    if (model.status != "partial") return;
    data.local = census::local_models(model);
    data.sign = census::functional_equation_sign(*data.local, model.q);
    data.torsion = census::complete_torsion(model, *data.local);
    data.lfun = counter.compute(model, *data.local, *data.sign);
    data.derived = census::derived_invariants(model, *data.local, *data.torsion, *data.lfun);
    data.candidates = census::admissible_candidates(*data.derived, data.local->fibres);
}

bool invariants(std::string_view backend, bool diagnostic) {
    census::LPolynomialCounter counter(backend);
    std::unique_ptr<census::CompressedRecords> records;
    bool success = true;
    if (!diagnostic) {
        records = std::make_unique<census::CompressedRecords>(std::cout);
        records->write(census::invariant_csv_header()+"\n");
    }
    if (diagnostic) std::cout << "q,key,kind,minimality,status,a1,a2,a3,a4,a6,"
                 "infinity_a1,infinity_a2,infinity_a3,infinity_a4,infinity_a6,"
                 "fibres,twist_fibres,degree,sign,height,twist_height,tamagawa_exponent,"
                 "torsion,torsion_points,lower_L,L_polynomial,trace_sums,"
                 "rank,geometric_rank,picard_number,geometric_picard_number,"
                 "P2,zeta_sign,point_counts,brauer_height,artin_tate_product,ns_splitting_degree,"
                 "j_degree,j_inseparability_exponent,sha_reg,"
                 "admissible_sha,admissible_regulators\n";
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(std::cin, line)) {
        ++line_number;
        try {
            std::istringstream input(line);
            std::string field, key, extra;
            if (!(input >> field)) continue;
            const auto q = unsigned_integer(field, "q");
            if (q != 2 && q != 3) throw std::invalid_argument("q must be the prime 2 or 3");
            if (!(input >> key) || input >> extra) throw std::invalid_argument("expected q key");
            const auto model = census::model_from_key(q, unsigned_integer(key, "key"));
            census::InvariantData data;
            std::optional<std::string> error;
            try {
                compute_invariants(model, counter, data);
            } catch (const std::exception& failure) {
                if (!records) throw;
                error = failure.what();
            }
            if (records) {
                records->write(census::invariant_csv(model, data, error)+"\n");
                if (error || (data.candidates && data.candidates->sha.empty())) success = false;
                continue;
            }
            std::cout << model.q << ',' << model.key << ',' << model.kind << ','
                      << model.minimality << ',' << model.status;
            write_coefficients(model.affine);
            write_coefficients(model.infinity);
            std::cout << ',';
            if (model.status == "partial") {
                const auto& local = *data.local;
                const auto& sign = *data.sign;
                const auto& torsion = *data.torsion;
                const auto& lfun = *data.lfun;
                census::write_fibres(std::cout, local.fibres, true);
                std::cout << ',';
                census::write_fibres(std::cout, local.twist_fibres, true);
                std::cout << ',' << sign.degree << ',' << sign.sign << ',' << sign.height
                          << ',' << sign.twist_height << ',' << sign.tamagawa_exponent;
                std::cout << ',';
                census::write_list(std::cout, torsion.factors);
                std::cout << ',';
                census::write_joined(std::cout, torsion.points, ';');
                std::cout << ',';
                census::write_list(std::cout, lfun.coefficients, 1, 1+sign.degree/2);
                std::cout << ',';
                census::write_list(std::cout, lfun.coefficients, 0, lfun.coefficients.size());
                std::cout << ',';
                census::write_list(std::cout, lfun.trace_sums, 0, lfun.trace_sums.size());
                std::cout << ',';
                census::write_derived_invariants(std::cout, *data.derived);
                std::cout << ',';
                census::write_candidate_invariants(std::cout, *data.candidates);
            } else std::cout << std::string(26, ',');
            std::cout << '\n';
        } catch (const std::exception& error) {
            if (records) records->finish();
            throw std::runtime_error("line " + std::to_string(line_number) + ": " + error.what());
        }
    }
    if (records) records->finish();
    if (std::cin.bad()) throw std::runtime_error("could not read standard input");
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("could not write standard output");
    return success;
}

void count() {
    std::map<std::pair<int, int>, census::TraceTable> tables;
    std::cout << "p,degree,parameter,trace,points\n";
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(std::cin, line)) {
        ++line_number;
        try {
            std::istringstream input(line);
            std::string prime, extension, parameter, extra;
            if (!(input >> prime)) continue;
            if (!(input >> extension >> parameter) || input >> extra) {
                throw std::invalid_argument("expected p degree parameter");
            }
            const auto p = unsigned_integer(prime, "p");
            const auto degree = unsigned_integer(extension, "degree");
            if ((p != 2 && p != 3) || degree < 1 || degree > 12) {
                throw std::invalid_argument("p must be 2 or 3 and degree must be between 1 and 12");
            }
            std::uint64_t size = 1;
            for (unsigned i = 0; i < degree; ++i) size *= p;
            const bool zero = parameter == "zero";
            const auto exponent = zero ? 0 : unsigned_integer(parameter, "parameter");
            if (exponent >= size - 1) throw std::invalid_argument("parameter exponent out of range");
            const auto& table = tables.try_emplace({int(p), int(degree)}, p, degree).first->second;
            const int trace = zero ? table.nodal : table.ordinary[exponent];
            std::cout << p << ',' << degree << ',' << parameter << ',' << trace << ','
                      << std::int64_t(size) + 1 - trace << '\n';
        } catch (const std::exception& error) {
            throw std::runtime_error("line " + std::to_string(line_number) + ": " + error.what());
        }
    }
    if (std::cin.bad()) throw std::runtime_error("could not read standard input");
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("could not write standard output");
}
}

int main(int argc, char* argv[]) {
    if (argc == 2) {
        const std::string_view argument(argv[1]);
        if (argument == "help" || argument == "--help" || argument == "-h") {
            print_help(std::cout);
            return 0;
        }
    }
    if (argc >= 2 && std::string_view(argv[1]) == "count") {
        if (argc != 2) {
            std::cerr << "count accepts no options\n";
            return 2;
        }
        try {
            count();
            return 0;
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    if (argc >= 2 && std::string_view(argv[1]) == "invariants") {
        try {
            std::string_view backend = "auto";
            bool diagnostic = false;
            for (int i = 2; i < argc; ++i) {
                const std::string_view option(argv[i]);
                if (option == "--diagnostic") diagnostic = true;
                else if (option == "--backend" && i+1 < argc) backend = argv[++i];
                else throw std::invalid_argument("invalid invariants option");
            }
            return invariants(backend, diagnostic) ? 0 : 1;
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    print_help(std::cerr);
    return 2;
}
