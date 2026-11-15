#include <NTL/ZZ.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <flint/flint.h>
#include <gmp.h>
#include <zstd.h>

#include <iostream>
#include <string_view>

namespace {
void print_help(std::ostream& output) {
    output << "Usage: census-native <command>\n\n"
              "Commands:\n"
              "  help       Show this help\n\n"
              "Options:\n"
              "  -h, --help  Show this help\n\n"
              "Arithmetic commands are not available in this build.\n";
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
    print_help(std::cerr);
    return 2;
}
