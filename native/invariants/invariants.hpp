#pragma once

#include "../models/model.hpp"

namespace census {

struct TorsionData {
    std::vector<int> factors;
    std::vector<std::string> points;
    std::vector<int> orders;
};

TorsionData complete_torsion(const Model& model, const LocalModels& local);

}
