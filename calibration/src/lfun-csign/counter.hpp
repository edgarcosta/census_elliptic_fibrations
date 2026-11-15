#ifndef LFUN_COUNTER_HPP
#define LFUN_COUNTER_HPP
#include "field.hpp"
#include <array>
using Poly=std::array<int,25>;
struct Surface {
    int p,D,sign;
    uint64_t key;
    std::array<Poly,5> a;
    Poly primary,delta,A,B;
    explicit Surface(int characteristic):p(characteristic){}
    void invariants();
};
void evaluate(const Field&,const Poly&,std::vector<uint32_t>&);
int64_t sum_ordinary3(const Field&,const std::vector<uint32_t>&,const std::vector<uint32_t>&,bool&);
#endif
