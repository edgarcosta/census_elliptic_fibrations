#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include "search.hpp"
namespace f2 {
using U = std::uint64_t;
int degree(U x) { return x ? 63 - __builtin_clzll(x) : -1; }
U rem(U a, U b) {
    int d = degree(b);
    while (a && degree(a) >= d) a ^= b << (degree(a) - d);
    return a;
}
U mul(U a, U b) {
    U r = 0;
    while (b) { if (b & 1) r ^= a; a <<= 1; b >>= 1; }
    return r;
}
U gcd(U a, U b) {
    while (b) { U r = rem(a,b); a=b; b=r; }
    return a;
}
U inverse(U a, U h) {
    U u=rem(a,h), v=h, s=1, t=0;
    if (!u) throw std::runtime_error("nonunit");
    while (u!=1) {
        int j=degree(u)-degree(v);
        if (j<0) { std::swap(u,v); std::swap(s,t); j=-j; }
        u ^= v << j; s ^= t << j;
    }
    return rem(s,h);
}
struct Place {
    U h; int d; std::vector<U> xs; std::array<std::array<unsigned char,256>,4> red;
    unsigned reduce(U f) const {
        return red[0][f&255] ^ red[1][(f>>8)&255] ^ red[2][(f>>16)&255] ^ red[3][(f>>24)&255];
    }
};
std::vector<Place> local_tables(const std::array<U,5>& aa) {
    std::vector<Place> ps;
    for(int d=1;d<=4;++d) for(U h=U(1)<<d;h<(U(1)<<(d+1));++h) {
        bool irr=true;
        for(int e=1;e<=d/2;++e) for(U f=U(1)<<e;f<(U(1)<<(e+1));++f)
            if(rem(h,f)==0) irr=false;
        if(!irr) continue;
        Place p; p.h=h; p.d=d;
        std::array<U,5> a; for(int i=0;i<5;++i) a[i]=rem(aa[i],h);
        for(U x=0;x<(U(1)<<d);++x) {
            U x2=rem(mul(x,x),h),L=rem(mul(a[0],x),h)^a[2];
            U R=rem(mul(x2,x)^mul(a[1],x2)^mul(a[3],x)^a[4],h);
            for(U y=0;y<(U(1)<<d);++y)
                if(rem(mul(y,y)^mul(L,y),h)==R) { p.xs.push_back(x); break; }
        }
        for(int j=0;j<4;++j) for(int c=0;c<256;++c) p.red[j][c]=rem(U(c)<<(8*j),h);
        ps.push_back(p);
    }
    return ps;
}
void emit(const std::array<U,5>& a,U A,U B,int b,std::ostream& out,bool& first,Stats& stats) {
    U B2=mul(B,B), B3=mul(B2,B);
    U R=mul(mul(A,A),A)^mul(a[1],mul(mul(A,A),B2))^mul(a[3],mul(A,mul(B2,B2)))^mul(a[4],mul(B3,B3));
    U L=mul(a[0],mul(A,B))^mul(a[2],B3);
    int m=7+3*b;
    std::array<U,64> basis{},sol{};
    U C=0;
    if(!L) {
        for(int j=0;j<64;++j) if((R>>j)&1) {
            if(j%2 || j/2>=m) return;
            C|=U(1)<<(j/2);
        }
    } else {
        for(int j=0;j<m;++j) {
            U v=(U(1)<<(2*j))^(L<<j), s=U(1)<<j;
            while(v) {
                int k=degree(v);
                if(!basis[k]) {basis[k]=v; sol[k]=s; break;}
                v^=basis[k];s^=sol[k];
            }
        }
        U v=R;
        while(v) {
            int k=degree(v); if(!basis[k]) return;
            v^=basis[k]; C^=sol[k];
        }
    }
    for(int j=0;j<(L?2:1);++j) {
        U y=j?C^L:C;
        if(!first) out<<',';
        first=false; out<<'['<<A<<','<<B<<','<<y<<']'; ++stats.points;
    }
}
void run(const std::array<U,5>& aa, int b, std::ostream& out, Stats& stats) {
    int N=2*b+4, target=N+1;
    auto places=local_tables(aa);
    U visited=0, passed=0;
    bool first=true;
    for(U g=1;g<(U(1)<<(b+1));++g) {
        U g2=mul(g,g); bool skip=false;
        std::vector<std::vector<U>> choices(places.size());
        std::vector<std::array<U,2>> masks(places.size());
        for(unsigned i=0;i<places.size();++i) {
            auto &p=places[i]; auto &xs=choices[i];
            if(rem(g,p.h)==0) {
                for(U x=1;x<(U(1)<<p.d);++x) xs.push_back(x);
            } else {
                for(U x:p.xs) xs.push_back(rem(mul(x,g2),p.h));
            }
            if(xs.empty()) { skip=true; break; }
            for(U x:xs) masks[i][x/64] |= U(1)<<(x%64);
        }
        if(skip) continue;
        const U inf=std::numeric_limits<U>::max()/2;
        std::vector<U> cost(target+1,inf);
        std::vector<std::vector<int>> paths(target+1); cost[0]=1;
        for(unsigned i=0;i<places.size();++i) {
            for(int n=target-places[i].d;n>=0;--n) {
                if(cost[n]==inf) continue;
                U c=cost[n]*choices[i].size(); int k=n+places[i].d;
                if(c<cost[k]) { cost[k]=c; paths[k]=paths[n]; paths[k].push_back(i); }
            }
        }
        if(cost[target]==inf) throw std::runtime_error("insufficient interpolation degree");
        U total=cost[target];
        auto selected=paths[target]; U M=1;
        for(int i:selected) M=mul(M,places[i].h);
        if(degree(M)!=target) throw std::runtime_error("CRT degree");
        std::vector<std::vector<U>> offsets;
        for(int i:selected) {
            U hi=places[i].h, cofactor=0, r=M;
            while(r && degree(r)>=degree(hi)) {
                int j=degree(r)-degree(hi); cofactor ^= U(1)<<j; r ^= hi<<j;
            }
            if(r) throw std::runtime_error("CRT division");
            U basis=mul(cofactor,inverse(cofactor,hi));
            std::vector<U> v;
            for(U x:choices[i]) v.push_back(rem(mul(basis,x),M));
            offsets.push_back(v);
        }
        std::vector<int> remaining;
        for(unsigned i=0;i<places.size();++i)
            if(std::find(selected.begin(),selected.end(),i)==selected.end()) remaining.push_back(i);
        std::stable_sort(remaining.begin(),remaining.end(),[&](int i,int j){
            return (choices[i].size()<<places[j].d)<(choices[j].size()<<places[i].d);
        });
        std::vector<unsigned> digits(offsets.size(),0); U f=0;
        for(auto &v:offsets) f ^= v[0];
        for(U k=0;k<total;++k) {
            ++visited; bool valid=true;
            for(int i:remaining) {
                unsigned residue=places[i].reduce(f);
                if(((masks[i][residue/64]>>(residue%64))&1)==0) { valid=false; break; }
            }
            if(valid && gcd(f,g)==1) { ++passed; emit(aa,f,g,b,out,first,stats); }
            for(unsigned j=0;j<offsets.size();++j) {
                f ^= offsets[j][digits[j]];
                if(++digits[j]<offsets[j].size()) { f ^= offsets[j][digits[j]]; break; }
                digits[j]=0; f ^= offsets[j][0];
            }
        }
    }
    stats.visited=visited; stats.candidates=passed;
}
}
void search2(const std::array<std::uint64_t,5>& a, int b, std::ostream& out, Stats& stats) { f2::run(a,b,out,stats); }
