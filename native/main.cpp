#include "lfun.hpp"
#include "model.hpp"
#include "invariants.hpp"

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
              "  --coefficients            Read q a1 a2 a3 a4 a6 as base-q codes\n"
              "  --chart 'a_i(t,1)'        Input chart for invariants (default a_i(t,1))\n\n"
              "Input: one q key pair per line, q=2 or 3, decimal uint64 key.\n"
              "Output: CSV with both charts, cell labels, local fibres, signed L, torsion and derived invariants.\n"
              "Coefficient input uses a2=0 over F_2 and a1=a3=0 over F_3.\n"
              "The key and sign are derived from the equation; no sign is supplied.\n"
              "Smooth everywhere minimal rows are partial; other cells are unsupported.\n";
    output << "\ncount [--table-cache DIRECTORY] reads p degree parameter rows.\n"
              "Parameter is an exponent 0 <= e < p^degree-1, or zero for parameter 0.\n"
              "The parameter at exponent e is alpha^e for the retained field polynomial.\n"
              "Degrees 1 through 12 are supported. Output: p,degree,parameter,trace,points.\n"
              "Equations: y^2+x*y=x^3+a over F_2; y^2=x^3+x^2+a over F_3.\n"
              "Cache files are checked against exact tables before counting.\n";
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
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            if (i != 0) std::cout << ' ';
            std::cout << coefficients[i];
        }
    }
}

void write_fibres(const std::vector<census::LocalFibre>& fibres) {
    if (fibres.empty()) std::cout << "[]";
    for (std::size_t i = 0; i < fibres.size(); ++i) {
        const auto& f = fibres[i];
        if (i != 0) std::cout << ';';
        std::cout << f.place << ':' << f.degree << ':' << f.type << ':' << f.index
                  << ':' << f.discriminant_valuation << ':' << f.conductor << ':'
                  << f.tamagawa << ':' << f.trace << ':';
        if (f.component_cycles.empty()) std::cout << "[]";
        for (std::size_t j = 0; j < f.component_cycles.size(); ++j) {
            if (j != 0) std::cout << ' ';
            std::cout << f.component_cycles[j];
        }
    }
}

template<class Values>
void write_list(const Values& values, std::size_t first, std::size_t end) {
    if (first == end) std::cout << "[]";
    for (std::size_t i = first; i < end; ++i) {
        if (i != first) std::cout << ' ';
        std::cout << values[i];
    }
}

void invariants(std::string_view backend, bool coefficients) {
    census::LPolynomialCounter counter(backend);
    std::cout << "q,key,kind,minimality,status,chart,a1,a2,a3,a4,a6,"
                 "infinity_a1,infinity_a2,infinity_a3,infinity_a4,infinity_a6,"
                 "fibres,twist_fibres,degree,sign,height,twist_height,tamagawa_exponent,"
                 "torsion,torsion_points,torsion_orders,lower_L,L_polynomial,trace_sums,"
                 "rank,geometric_rank,picard_number,geometric_picard_number,"
                 "P2,zeta_sign,point_counts,brauer_height,artin_tate_product,ns_splitting_degree,"
                 "j_degree,j_inseparability_exponent,is_isotrivial,sha_reg\n";
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
            std::uint64_t packed = 0;
            if (coefficients) {
                std::array<std::string, 5> codes;
                if (!(input >> codes[0] >> codes[1] >> codes[2] >> codes[3] >> codes[4])
                    || input >> extra) throw std::invalid_argument("expected q a1 a2 a3 a4 a6");
                for (int j = 0; j < 5; ++j) {
                    const auto code = unsigned_integer(codes[j], "coefficient code");
                    const bool included = q == 2 ? j != 1 : j != 0 && j != 2;
                    if (!included && code != 0) {
                        throw std::invalid_argument("coefficients must use the census normal form");
                    }
                    std::uint64_t bound = 1;
                    for (int k = 0; k <= census::model_degrees[j]; ++k) bound *= q;
                    if (code >= bound) throw std::invalid_argument("coefficient exceeds K3 degree bound");
                    if (included) packed = packed*bound+code;
                }
            } else {
                if (!(input >> key) || input >> extra) throw std::invalid_argument("expected q key");
                packed = unsigned_integer(key, "key");
            }
            const auto model = census::model_from_key(q, packed);
            census::LocalModels local;
            census::SignData sign{};
            census::TorsionData torsion;
            census::LPolynomial lfun;
            census::DerivedInvariants derived;
            if (model.status == "partial") {
                local = census::local_models(model);
                sign = census::functional_equation_sign(local, model.q);
                torsion = census::complete_torsion(model, local);
                lfun = counter.compute(model, local, sign);
                derived = census::derived_invariants(model, local, torsion, lfun);
            }
            std::cout << model.q << ',' << model.key << ',' << model.kind << ','
                      << model.minimality << ',' << model.status << ",\""
                      << census::model_chart << '\"';
            write_coefficients(model.affine);
            write_coefficients(model.infinity);
            std::cout << ',';
            if (model.status == "partial") {
                write_fibres(local.fibres);
                std::cout << ',';
                write_fibres(local.twist_fibres);
                std::cout << ',' << sign.degree << ',' << sign.sign << ',' << sign.height
                          << ',' << sign.twist_height << ',' << sign.tamagawa_exponent;
                std::cout << ',';
                if (torsion.factors.empty()) std::cout << "[]";
                for (std::size_t i = 0; i < torsion.factors.size(); ++i) {
                    if (i) std::cout << ' ';
                    std::cout << torsion.factors[i];
                }
                std::cout << ',';
                for (std::size_t i = 0; i < torsion.points.size(); ++i) {
                    if (i) std::cout << ';';
                    std::cout << torsion.points[i];
                }
                std::cout << ',';
                for (std::size_t i = 0; i < torsion.orders.size(); ++i) {
                    if (i) std::cout << ' ';
                    std::cout << torsion.orders[i];
                }
                std::cout << ',';
                write_list(lfun.coefficients, 1, 1+sign.degree/2);
                std::cout << ',';
                write_list(lfun.coefficients, 0, lfun.coefficients.size());
                std::cout << ',';
                write_list(lfun.trace_sums, 0, lfun.trace_sums.size());
                std::cout << ',' << derived.rank << ',' << derived.geometric_rank
                          << ',' << derived.picard_number << ',' << derived.geometric_picard_number << ',';
                write_list(derived.p2, 0, derived.p2.size());
                std::cout << ',' << derived.zeta_sign << ',';
                write_list(derived.point_counts, 0, derived.point_counts.size());
                std::cout << ',';
                if (derived.brauer_height) std::cout << derived.brauer_height;
                else std::cout << "infinity";
                std::cout << ',' << derived.artin_tate_product << ',' << derived.ns_splitting_degree
                          << ',' << derived.j_degree << ',' << derived.j_inseparability_exponent
                          << ',' << (derived.j_degree == 0 ? "true" : "false") << ',';
                std::cout << numerator(derived.sha_reg);
                if (denominator(derived.sha_reg) != 1) std::cout << '/' << denominator(derived.sha_reg);
            } else std::cout << std::string(26, ',');
            std::cout << '\n';
        } catch (const std::exception& error) {
            throw std::runtime_error("line " + std::to_string(line_number) + ": " + error.what());
        }
    }
    if (std::cin.bad()) throw std::runtime_error("could not read standard input");
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("could not write standard output");
}

void count(const std::filesystem::path& cache) {
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
            const auto& table = tables.try_emplace({int(p), int(degree)}, p, degree, cache).first->second;
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
        if (argc != 2 && !(argc == 4 && std::string_view(argv[2]) == "--table-cache"
                          && std::string_view(argv[3]) != "")) {
            std::cerr << "count accepts --table-cache DIRECTORY\n";
            return 2;
        }
        try {
            count(argc == 4 ? argv[3] : "");
            return 0;
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    if (argc >= 2 && std::string_view(argv[1]) == "invariants") {
        try {
            std::string_view backend = "auto";
            bool coefficients = false;
            for (int i = 2; i < argc; ++i) {
                const std::string_view option(argv[i]);
                if (option == "--coefficients") coefficients = true;
                else if (option == "--backend" && i+1 < argc) backend = argv[++i];
                else if (option == "--chart" && i+1 < argc
                         && std::string_view(argv[i+1]) == census::model_chart) ++i;
                else throw std::invalid_argument("invariants requires chart 'a_i(t,1)' and valid options");
            }
            invariants(backend, coefficients);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    print_help(std::cerr);
    return 2;
}
