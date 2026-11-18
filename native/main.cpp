#include "models/model.hpp"
#include "lfunction/lfun.hpp"

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
              "  invariants  Decode and classify q/key rows from standard input\n"
              "  count       Count points on the ordinary cubic families\n\n"
              "Options:\n"
              "  -h, --help                Show this help\n"
              "  --diagnostic              Emit model comparison CSV (default)\n"
              "Input: one q key pair per line, q=2 or 3, decimal uint64 key.\n"
              "Output: CSV with both charts, cell labels and arithmetic status.\n"
              "Smooth everywhere minimal rows are partial; other cells are unsupported.\n";
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
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            if (i != 0) std::cout << ' ';
            std::cout << coefficients[i];
        }
    }
}

void invariants() {
    std::cout << "q,key,kind,minimality,status,a1,a2,a3,a4,a6,"
                 "infinity_a1,infinity_a2,infinity_a3,infinity_a4,infinity_a6\n";
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(std::cin, line)) {
        ++line_number;
        try {
            std::istringstream input(line);
            std::string field, key, extra;
            if (!(input >> field)) continue;
            if (!(input >> key) || input >> extra) {
                throw std::invalid_argument("expected q key");
            }
            const auto q = unsigned_integer(field, "q");
            if (q != 2 && q != 3) throw std::invalid_argument("q must be the prime 2 or 3");
            const auto model = census::model_from_key(q, unsigned_integer(key, "key"));
            std::cout << model.q << ',' << model.key << ',' << model.kind << ','
                      << model.minimality << ',' << model.status;
            write_coefficients(model.affine);
            write_coefficients(model.infinity);
            std::cout << '\n';
        } catch (const std::exception& error) {
            throw std::runtime_error("line " + std::to_string(line_number) + ": " + error.what());
        }
    }
    if (std::cin.bad()) throw std::runtime_error("could not read standard input");
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("could not write standard output");
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
            for (int i = 2; i < argc; ++i) {
                const std::string_view option(argv[i]);
                if (option == "--diagnostic") continue;
                throw std::invalid_argument("invalid invariants option");
            }
            invariants();
            return 0;
        } catch (const std::exception& error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    print_help(std::cerr);
    return 2;
}
