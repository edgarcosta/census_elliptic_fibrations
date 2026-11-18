#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace census {

inline constexpr std::array<int, 5> model_degrees = {2, 4, 6, 8, 12};
inline constexpr const char* model_chart = "a_i(t,1)";

struct Model {
    unsigned q;
    std::uint64_t key;
    std::array<std::vector<int>, 5> affine;
    std::array<std::vector<int>, 5> infinity;
    std::string kind;
    std::string minimality;
    std::string status;
};

struct LocalFibre {
    std::string place;
    int degree;
    std::string type;
    int index;
    int discriminant_valuation;
    int conductor;
    int tamagawa;
    int trace;
    std::vector<int> component_cycles;
};

struct LocalModels {
    std::vector<LocalFibre> fibres;
    std::vector<LocalFibre> twist_fibres;
};

Model model_from_key(unsigned q, std::uint64_t key);
LocalModels local_models(const Model& model);

}
