#pragma once

#include "field.hpp"
#include "model.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string_view>
#include <vector>

namespace census {

struct TraceTable {
    TraceField field;
    std::vector<int> ordinary;
    int nodal;

    TraceTable(int p, int degree, const std::filesystem::path& cache = {});
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
    std::vector<boost::multiprecision::cpp_int> coefficients;
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
