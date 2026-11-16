#include "invariants.hpp"

#include <iostream>

// Prime-field keys cannot express the constant-field controls in the Magma API.
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
    census::write_derived_invariants(std::cout, census::derived_invariants(model, local, {}, L, 12));
    std::cout << '\n';
    model.affine[3].push_back(1);
    model.affine[3][6] = 0;
    local.fibres[0].tamagawa = 2;
    local.fibres[0].component_cycles = {1};
    // LFunction after squaring Frobenius on y^2+y=x^3+t^7*x.
    L.coefficients = {1, -4, 0, 240, -960, 0, 24576, -98304, 0, 1703936,
                      -6815744, 0, 100663296, -402653184, 0, 4026531840LL,
                      -16106127360LL, 0, 68719476736LL, -274877906944LL};
    census::write_derived_invariants(std::cout, census::derived_invariants(model, local, {}, L, 2));
    std::cout << '\n';
}
