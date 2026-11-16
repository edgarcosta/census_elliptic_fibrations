#include "invariants.hpp"

#include <iostream>

// Prime-field keys cannot express the constant-field controls in the Magma API.
void write(const census::DerivedInvariants& r) {
    const auto list = [](const auto& values) {
        for (std::size_t i = 0; i < values.size(); ++i) std::cout << (i ? " " : "") << values[i];
    };
    std::cout << r.rank << ',' << r.geometric_rank << ',' << r.picard_number
              << ',' << r.geometric_picard_number << ',';
    list(r.p2);
    std::cout << ',' << r.zeta_sign << ',';
    list(r.point_counts);
    std::cout << ',';
    if (r.brauer_height) std::cout << r.brauer_height;
    else std::cout << "infinity";
    std::cout << ',' << r.artin_tate_product << ',' << r.ns_splitting_degree
              << ',' << r.j_degree << ',' << r.j_inseparability_exponent << ','
              << (r.j_degree == 0 ? "true" : "false") << ',' << numerator(r.sha_reg);
    if (denominator(r.sha_reg) != 1) std::cout << '/' << denominator(r.sha_reg);
    std::cout << '\n';
}

int main() {
    census::Model model{};
    model.q = 2;
    model.status = "partial";
    model.affine = {{{0}, {0}, {1}, {0, 0, 0, 0, 0, 0, 1}, {0}}};
    census::LocalModels local;
    census::LocalFibre fibre{};
    fibre.tamagawa = 4;
    fibre.component_cycles.assign(13, 1);
    local.fibres.push_back(fibre);
    census::LPolynomial L;
    L.coefficients = {1};
    for (int n = 0; n < 7; ++n) {
        const int degree = L.coefficients.size()-1;
        L.coefficients.push_back(0);
        for (int i = degree; i >= 0; --i) L.coefficients[i+1] -= 4096*L.coefficients[i];
    }
    write(census::derived_invariants(model, local, {}, L, 12));
    model.affine[3].push_back(1);
    model.affine[3][6] = 0;
    local.fibres[0].tamagawa = 2;
    local.fibres[0].component_cycles = {1};
    // LFunction after squaring Frobenius on y^2+y=x^3+t^7.
    L.coefficients = {1, -4, 0, 240, -960, 0, 24576, -98304, 0, 1703936,
                      -6815744, 0, 100663296, -402653184, 0, 4026531840LL,
                      -16106127360LL, 0, 68719476736LL, -274877906944LL};
    write(census::derived_invariants(model, local, {}, L, 2));
}
