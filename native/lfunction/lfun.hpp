#pragma once

#include "field.hpp"
#include "../models/model.hpp"

#include <gmpxx.h>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace census {

struct TraceTable {
    TraceField field;
    std::vector<int> ordinary;
    int nodal;

    TraceTable(int p, int degree);
};

struct SignData {
    int degree;
    int sign;
    int height;
    int twist_height;
    int tamagawa_exponent;
};

SignData functional_equation_sign(const LocalModels& local, unsigned q);

struct LPolynomial {
    std::vector<mpz_class> coefficients;
    std::vector<std::int64_t> trace_sums;
};

class LPolynomialCounter {
public:
    explicit LPolynomialCounter(std::string_view backend = "auto");
    ~LPolynomialCounter();
    LPolynomial compute(const Model& model, const LocalModels& local, const SignData& sign);

private:
    struct State;
    std::unique_ptr<State> state;
};

}
