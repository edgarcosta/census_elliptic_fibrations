#pragma once

#include "field.hpp"
#include "model.hpp"

#include <filesystem>
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

}
