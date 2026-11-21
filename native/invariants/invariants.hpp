#pragma once

#include "../models/model.hpp"

namespace census {

struct TorsionData {
    std::vector<int> factors;
    std::vector<std::string> points;
};

TorsionData complete_torsion(const Model& model, const LocalModels& local);

}
