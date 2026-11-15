#ifndef LFUN_TATE_HPP
#define LFUN_TATE_HPP
#include "counter.hpp"
#include <string>
#include <vector>
struct LocalData {
    std::string place,kodaira;
    int degree,vd,conductor,tamagawa,trace,scalings;
};
struct SignData {
    int sign=0,height=0,twist_height=0,tamagawa_exponent=0;
};
std::vector<LocalData> tate(const Surface&,bool* everywhere_minimal=nullptr,SignData* sign=nullptr);
#endif
