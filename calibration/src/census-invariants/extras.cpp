#include <algorithm>
#include <numeric>
#include <NTL/lzz_pX.h>
#include <array>
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace NTL;
double cpu() { timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t); return t.tv_sec+1e-9*t.tv_nsec; }
struct Row { long long key; zz_pX c4,delta; int degree=0,insep=0; };
zz_pX decode(long long n,int q) { zz_pX a; int i=0; while(n) { SetCoeff(a,i++,n%q); n/=q; } return a; }
void jdata(Row& r,int q) {
    zz_pX n=power(r.c4,3),d=r.delta,g=GCD(n,d);
    n/=g; d/=g; r.degree=std::max(deg(n),deg(d)); r.insep=0;
    if(r.degree>0) {
        int common=0;
        for(int i=1;i<=deg(n);++i) if(!IsZero(coeff(n,i))) common=std::gcd(common,i);
        for(int i=1;i<=deg(d);++i) if(!IsZero(coeff(d,i))) common=std::gcd(common,i);
        while(common%q==0) { ++r.insep; common/=q; }
    }
}
int quasi(long long key) {
    std::array<int,13> b{}; bool cube=true;
    for(int i=0;i<13;++i) { b[i]=key%3; key/=3; if(i%3 && b[i]) cube=false; }
    if(cube) return 0;
    auto nonminimal=[](const std::array<int,13>& a) { return !(a[1]||a[2]||a[4]||a[5]); };
    auto reverse=b; std::reverse(reverse.begin(),reverse.end());
    if(nonminimal(reverse)) return 1;
    for(int v=0;v<3;++v) {
        std::array<int,13> a{};
        for(int i=12;i>=0;--i) {
            for(int j=12;j>=1;--j) a[j]=(a[j]*v+a[j-1])%3;
            a[0]=(a[0]*v+b[i])%3;
        }
        if(nonminimal(a)) return 1;
    }
    return 2;
}
volatile long long sink=0;
int main(int argc,char**argv) {
    assert(argc==3); int q=std::stoi(argv[1]); zz_p::init(q);
    std::ifstream f(argv[2]); std::vector<Row> rows;
    int p,sg; long long key; std::array<long long,5> code;
    while(f>>p>>key>>sg) {
        assert(p==q); std::array<zz_pX,5> a;
        for(int i=0;i<5;++i) { f>>code[i]; a[i]=decode(code[i],q); }
        auto b2=a[0]*a[0]+4*a[1], b4=a[0]*a[2]+2*a[3], b6=a[2]*a[2]+4*a[4];
        auto b8=a[0]*a[0]*a[4]+4*a[1]*a[4]-a[0]*a[2]*a[3]+a[1]*a[2]*a[2]-a[3]*a[3];
        Row r; r.key=key; r.c4=b2*b2-24*b4;
        r.delta=-b2*b2*b8-8*b4*b4*b4-27*b6*b6+9*b2*b4*b6;
        jdata(r,q); rows.push_back(r);
    }
    for(auto r:rows) std::cout<<r.key<<","<<r.degree<<","<<r.insep<<"\n";
}
