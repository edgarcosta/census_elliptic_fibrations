#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include "search.hpp"
namespace f3 {
using U=std::uint64_t;
// Disjoint masks record coefficients 1 and 2; addition has no carries between powers.
struct P { U one=0,two=0; };
P add(P a,P b) { U za=~(a.one|a.two),zb=~(b.one|b.two); return {(za&b.one)|(a.one&zb)|(a.two&b.two),(za&b.two)|(a.two&zb)|(a.one&b.one)}; }
P neg(P a) { return {a.two,a.one}; }
P shift(P a,int j) { return {a.one<<j,a.two<<j}; }
int degree(P a) { U bits=a.one|a.two; return bits?63-__builtin_clzll(bits):-1; }
int leading(P a) { return (a.two>>degree(a))&1?2:1; }
P fromcode(U n) { P a; int i=0; while(n) { int c=n%3; if(c==1)a.one|=U(1)<<i; if(c==2)a.two|=U(1)<<i; ++i;n/=3; }return a; }
U code(P a) { U n=0,p=1;while(a.one|a.two){n+=p*((a.one&1)+2*(a.two&1));p*=3;a.one>>=1;a.two>>=1;}return n; }
P mul(P a,P b) { P r;while(b.one|b.two){if(b.one&1)r=add(r,a);if(b.two&1)r=add(r,neg(a));a=shift(a,1);b.one>>=1;b.two>>=1;}return r; }
std::pair<P,P> divide(P a,P b) {
    P q; int d=degree(b); if(d<0)throw std::runtime_error("division by zero");
    while(degree(a)>=d) { int j=degree(a)-d,c=leading(a)==leading(b)?1:2; P v=shift(b,j);a=add(a,c==1?neg(v):v);if(c==1)q.one|=U(1)<<j;else q.two|=U(1)<<j; }
    return {q,a};
}
P rem(P a,P b) { return divide(a,b).second; }
P gcd(P a,P b) { while(degree(b)>=0){P r=rem(a,b);a=b;b=r;}return a; }
P inverse(P a,P h) {
    P r0=h,r1=rem(a,h),s0,s1{1,0};
    if(degree(r1)<0)throw std::runtime_error("nonunit");
    while(degree(r1)>0) { auto qr=divide(r0,r1);P s=add(s0,neg(mul(qr.first,s1)));r0=r1;r1=qr.second;s0=s1;s1=s; }
    return rem(leading(r1)==1?s1:neg(s1),h);
}
struct Place { P h;int d,q;std::vector<int>xs;std::array<std::array<unsigned char,256>,8> red; };
std::vector<Place> local_tables(const std::array<U,5>& aa,int chunks) {
    std::vector<Place> ps;
    U pd=1;
    for(int d=1;d<=4;++d) {
        pd*=3;
        for(U n=0;n<pd;++n) {
            P h=fromcode(pd+n); bool irr=true; U pe=1;
            for(int e=1;e<=d/2;++e) {pe*=3;for(U f=pe;f<2*pe;++f) if(degree(rem(h,fromcode(f)))<0) irr=false;}
            if(!irr) continue;
            Place p;p.h=h;p.d=d;p.q=pd;
            std::array<P,5> a;for(int i=0;i<5;++i)a[i]=rem(fromcode(aa[i]),h);
            std::vector<bool> squares(pd,false);
            for(U y=0;y<pd;++y)squares[code(rem(mul(fromcode(y),fromcode(y)),h))]=true;
            for(U x=0;x<pd;++x) {
                P X=fromcode(x),x2=rem(mul(X,X),h),L=add(rem(mul(a[0],X),h),a[2]);
                P R=add(add(mul(x2,X),mul(a[1],x2)),add(mul(a[3],X),a[4]));
                if(squares[code(rem(add(mul(L,L),R),h))])p.xs.push_back(x);
            }
            for(int j=0;j<chunks;++j)for(int n=0;n<256;++n){U a=n&15,c=n>>4;p.red[j][n]=(a&c)?0:code(rem(shift(P{a,c},4*j),p.h));}
            ps.push_back(p);
        }
    }
    return ps;
}
void emit(const std::array<U,5>& aa,P A,P B,int b,std::ostream& out,bool& first,Stats& stats) {
    std::array<P,5>a;for(int i=0;i<5;++i)a[i]=fromcode(aa[i]);
    P B2=mul(B,B),B3=mul(B2,B),A2=mul(A,A);
    P R=add(add(mul(A2,A),mul(a[1],mul(A2,B2))),add(mul(a[3],mul(A,mul(B2,B2))),mul(a[4],mul(B3,B3))));
    P L=add(mul(a[0],mul(A,B)),mul(a[2],B3));
    P H=add(mul(L,L),R),C;
    int d=degree(H);
    if(d>=0) {
        if(d%2 || leading(H)!=1) return;
        int n=d/2; C.one=U(1)<<n;
        for(int k=n-1;k>=0;--k) {
            P diff=add(H,neg(mul(C,C)));
            int j=n+k, c=(diff.one>>j&1)+2*(diff.two>>j&1);
            if(c==1) C.two|=U(1)<<k;
            if(c==2) C.one|=U(1)<<k;
        }
        if(degree(add(mul(C,C),neg(H)))>=0) return;
    }
    for(int j=0;j<(d>=0?2:1);++j) {
        P y=add(L,j?C:neg(C));
        if(degree(y)>6+3*b) throw std::runtime_error("root degree");
        if(!first)out<<',';
        first=false;out<<'['<<code(A)<<','<<code(B)<<','<<code(y)<<']';++stats.points;
    }
}
void run(const std::array<U,5>& aa, int b, std::ostream& out, Stats& stats) {
    int N=2*b+4,target=N+1,chunks=(N+4)/4;
    std::array<U,32> powers; powers[0]=1;for(int i=1;i<32;++i)powers[i]=3*powers[i-1];
    std::array<std::vector<unsigned char>,6> addition;
    for(int d=1;d<=4;++d){int q=powers[d];addition[d].resize(q*q);for(int x=0;x<q;++x)for(int y=0;y<q;++y)addition[d][x*q+y]=code(add(fromcode(x),fromcode(y)));}
    auto places=local_tables(aa,chunks);
    U visited=0,passed=0;bool first=true;
    for(int dg=0;dg<=b;++dg)for(U ng=0;ng<powers[dg];++ng){
        U gcode=powers[dg]+ng;P g=fromcode(gcode),g2=mul(g,g);bool skip=false;
        std::vector<std::vector<int>>choices(places.size());std::vector<std::array<U,4>>masks(places.size());
        for(unsigned i=0;i<places.size();++i){auto &p=places[i];auto &xs=choices[i];
            if(degree(rem(g,p.h))<0){for(int x=1;x<p.q;++x)xs.push_back(x);}else{for(int x:p.xs)xs.push_back(code(rem(mul(fromcode(x),g2),p.h)));}
            if(xs.empty()){skip=true;break;}for(int x:xs)masks[i][x/64]|=U(1)<<(x%64);
        }
        if(skip) continue;
        const U inf=std::numeric_limits<U>::max()/2;std::vector<U>cost(target+1,inf);std::vector<std::vector<int>>paths(target+1);cost[0]=1;
        for(unsigned i=0;i<places.size();++i){for(int n=target-places[i].d;n>=0;--n){if(cost[n]==inf)continue;U c=cost[n]*choices[i].size();int k=n+places[i].d;if(c<cost[k]){cost[k]=c;paths[k]=paths[n];paths[k].push_back(i);}}}
        if(cost[target]==inf)throw std::runtime_error("insufficient interpolation degree");
        auto selected=paths[target];P M{1,0};for(int i:selected)M=mul(M,places[i].h);
        if(degree(M)!=target) throw std::runtime_error("CRT degree");
        std::vector<std::vector<P>> offsets;
        for(int i:selected){P hi=places[i].h;auto qr=divide(M,hi);if(degree(qr.second)>=0)throw std::runtime_error("CRT division");P basis=mul(qr.first,inverse(qr.first,hi));std::vector<P>v;for(int x:choices[i])v.push_back(rem(mul(basis,fromcode(x)),M));offsets.push_back(v);}
        std::vector<int>remaining;for(unsigned i=0;i<places.size();++i)if(std::find(selected.begin(),selected.end(),i)==selected.end())remaining.push_back(i);
        std::stable_sort(remaining.begin(),remaining.end(),[&](int i,int j){return choices[i].size()*places[j].q<choices[j].size()*places[i].q;});
        std::vector<unsigned>digits(offsets.size(),0);P f;for(auto &v:offsets)f=add(f,v[0]);
        for(U k=0;k<cost[target];++k){
            ++visited;bool valid=true;
            for(int i:remaining){auto &p=places[i];int r=0;for(int j=0;j<chunks;++j){int n=((f.one>>(4*j))&15)|(((f.two>>(4*j))&15)<<4);r=addition[p.d][r*p.q+p.red[j][n]];}if(((masks[i][r/64]>>(r%64))&1)==0){valid=false;break;}}
            if(valid && degree(gcd(f,g))==0){++passed;emit(aa,f,g,b,out,first,stats);}
            for(unsigned j=0;j<offsets.size();++j){f=add(f,neg(offsets[j][digits[j]]));if(++digits[j]<offsets[j].size()){f=add(f,offsets[j][digits[j]]);break;}digits[j]=0;f=add(f,offsets[j][0]);}
        }
    }
    stats.visited=visited; stats.candidates=passed;
}
}
void search3(const std::array<std::uint64_t,5>& a, int b, std::ostream& out, Stats& stats) { f3::run(a,b,out,stats); }
