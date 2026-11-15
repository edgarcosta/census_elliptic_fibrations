#include "model.hpp"

#include <charconv>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
void print_help(std::ostream& output) {
    output << "Usage: census-native <command>\n\n"
              "Commands:\n"
              "  help        Show this help\n"
              "  invariants  Decode and classify q/key rows from standard input\n\n"
              "Options:\n"
              "  -h, --help                Show this help\n"
              "  --chart 'a_i(t,1)'        Input chart for invariants (default a_i(t,1))\n\n"
              "Input: one q key pair per line, q=2 or 3, decimal uint64 key.\n"
              "Output: CSV with both charts, cell labels and arithmetic status.\n"
              "Smooth everywhere minimal rows are partial; other cells are unsupported.\n";
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
    std::cout << "q,key,kind,minimality,status,chart,a1,a2,a3,a4,a6,"
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
                      << model.minimality << ',' << model.status << ",\""
                      << census::model_chart << '\"';
            write_coefficients(model.affine);
            write_coefficients(model.infinity);
            std::cout << '\n';
        } catch (const std::invalid_argument& error) {
            throw std::invalid_argument("line " + std::to_string(line_number) + ": " + error.what());
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
    if (argc >= 2 && std::string_view(argv[1]) == "invariants") {
        if (argc != 2 && !(argc == 4 && std::string_view(argv[2]) == "--chart"
                          && std::string_view(argv[3]) == census::model_chart)) {
            std::cerr << "invariants requires chart 'a_i(t,1)'\n";
            return 2;
        }
        try {
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
