#include "../../native/invariants/invariants.hpp"

#include <iostream>
#include <sstream>

int main() {
    std::cout << "admissible_sha,admissible_regulators,split_method,regulator,sha,brauer_order,ns_discriminant\n";
    std::string line;
    std::getline(std::cin, line);
    while (std::getline(std::cin, line)) {
        std::istringstream input(line);
        std::string rank, rho, n, product, rows;
        std::getline(input, rank, ',');
        std::getline(input, rho, ',');
        std::getline(input, n, ',');
        std::getline(input, product, ',');
        std::getline(input, rows, ',');
        census::DerivedInvariants data{};
        data.rank = std::stoi(rank);
        data.picard_number = std::stoi(rho);
        data.artin_tate_product = mpz_class(n);
        data.sha_reg = mpq_class(product);
        data.sha_reg.canonicalize();
        std::vector<census::LocalFibre> fibres;
        if (rows != "[]") {
            std::istringstream entries(rows);
            std::string entry;
            while (std::getline(entries, entry, ';')) {
                for (auto& c : entry) if (c == ':') c = ' ';
                std::istringstream fields(entry);
                census::LocalFibre f{};
                fields >> f.type >> f.index >> f.tamagawa >> f.degree;
                fibres.push_back(f);
            }
        }
        census::write_candidate_invariants(std::cout, data, census::admissible_candidates(data, fibres));
        std::cout << '\n';
    }
}
