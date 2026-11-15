#pragma once
#include <array>
#include <cstdint>
#include <ostream>
struct Stats { std::uint64_t visited=0,candidates=0,points=0; };
void search2(const std::array<std::uint64_t,5>&,int,std::ostream&,Stats&);
void search3(const std::array<std::uint64_t,5>&,int,std::ostream&,Stats&);
