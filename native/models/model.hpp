#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace census {

inline constexpr std::array<int, 5> model_degrees = {2, 4, 6, 8, 12};

struct Model {
    unsigned q;
    std::uint64_t key;
    std::array<std::vector<int>, 5> affine;
    std::array<std::vector<int>, 5> infinity;
    std::string kind;
    std::string minimality;
    std::string status;
};

Model model_from_key(unsigned q, std::uint64_t key);

}
