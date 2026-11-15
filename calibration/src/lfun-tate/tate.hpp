#ifndef LFUN_TATE_HPP
#define LFUN_TATE_HPP
#include "counter.hpp"
#include <string>
#include <vector>
struct LocalData {
    std::string place,kodaira;
    int degree,vd,conductor,tamagawa,trace,scalings;
};
std::vector<LocalData> tate(const Surface&,bool* everywhere_minimal=nullptr);
#endif
