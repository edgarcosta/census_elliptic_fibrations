#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
using I = int64_t;
using W = __int128_t;
using Poly = std::vector<I>;
#include "cyclotomics.inc"

I narrow(W n) { assert(n >= std::numeric_limits<I>::min() && n <= std::numeric_limits<I>::max()); return I(n); }
I power(I q,int n) { I a=1; while(n--) a*=q; return a; }
struct Rat {
    I n=0,d=1;
    Rat(W a=0,W b=1) {
        if(b<0) { a=-a; b=-b; }
        W x=a<0?-a:a,y=b;
        while(y) { W z=x%y; x=y; y=z; }
        n=narrow(a/x); d=narrow(b/x);
    }
    std::string str() const { return std::to_string(n)+(d==1?"":"/"+std::to_string(d)); }
};
struct Fiber { int d,type,n,c; };
struct Row {
    int q,timed,known; I key,tor,sha;
    Poly L,P2,quotient; std::vector<Fiber> fibers;
    std::vector<int> cycles;
    int rank=0,geom=0,rho=0,rhog=0,height=0,sign=0,tor_trivial=0;
    I nsdegree=1; Rat brdisc,sha_reg,disc,reg;
    std::array<I,5> points{};
};

std::vector<int> local_cycles(const Fiber& f) {
    switch(f.type) {
    case 1: return {};
    case 2: return {1};
    case 3: return f.c==3 ? std::vector<int>{1,1}:std::vector<int>{2};
    case 4: return f.c==3 ? std::vector<int>(6,1):std::vector<int>{1,1,2,2};
    case 5: return std::vector<int>(7,1);
    case 6: return std::vector<int>(8,1);
    case 7:
        if(f.n<=1) return {};
        if(f.c==f.n) return std::vector<int>(f.n-1,1);
        { std::vector<int> a((f.n-1)/2,2); if(f.n%2==0) a.push_back(1); return a; }
    case 8:
        if(f.n==0 && f.c==1) return {1,3};
        if(f.c==4) return std::vector<int>(f.n+4,1);
        { std::vector<int> a(f.n+2,1); a.push_back(2); return a; }
    }
    assert(false); return {};
}
void components(Row& r) {
    r.cycles.clear();
    for(auto f:r.fibers) for(int k:local_cycles(f)) r.cycles.push_back(k*f.d);
    std::sort(r.cycles.begin(),r.cycles.end());
}
void mult_binomial(Poly& a,int k,int q) {
    int d=a.size()-1; I b=power(q,k); a.resize(d+k+1);
    for(int i=d;i>=0;--i) a[i+k]=narrow(W(a[i+k])-W(b)*a[i]);
}
void p2(Row& r) {
    r.P2=r.L; mult_binomial(r.P2,1,r.q); mult_binomial(r.P2,1,r.q);
    for(int k:r.cycles) mult_binomial(r.P2,k,r.q);
    assert(r.P2.size()==23);
}
bool divide(const Poly& a,const Poly& b,Poly& quo) {
    if(a.size()<b.size()) return false;
    int n=a.size(),d=b.size()-1;
    std::array<W,23> work{}; for(int i=0;i<n;++i) work[i]=a[i];
    std::array<W,23> qq{};
    for(int i=0;i<n-d;++i) {
        qq[i]=work[i]/b[0];
        for(int j=0;j<=d;++j) work[i+j]-=qq[i]*b[j];
    }
    for(int i=n-d;i<n;++i) if(work[i]) return false;
    quo.resize(n-d); for(int i=0;i<n-d;++i) quo[i]=narrow(qq[i]);
    return true;
}
void rational_rank(Row& r) {
    r.quotient=r.L; r.rank=0; Poly z;
    while(divide(r.quotient,{1,-r.q},z)) { r.quotient=z; ++r.rank; }
    r.rho=2+r.rank+r.cycles.size();
}
std::array<std::vector<std::pair<int,Poly>>,4> scaled;
void geometric(Row& r) {
    Poly rem=r.L,z; r.geom=0; r.nsdegree=1;
    for(const auto& entry:scaled[r.q]) {
        const auto& c=entry.second;
        if(c.size()>rem.size()) continue;
        while(divide(rem,c,z)) {
            rem=z; r.geom+=c.size()-1;
            r.nsdegree=std::lcm(r.nsdegree,I(entry.first));
        }
    }
    r.rhog=2+r.geom;
    for(int k:r.cycles) { r.rhog+=k; r.nsdegree=std::lcm(r.nsdegree,I(k)); }
}
void height(Row& r) {
    int num=1,den=1;
    for(int i=1;i<int(r.L.size());++i) if(r.L[i]) {
        I a=r.L[i]<0?-r.L[i]:r.L[i]; int v=0;
        while(a%r.q==0) { a/=r.q; ++v; }
        if(v*den<num*i) { num=v; den=i; }
    }
    r.height=num==den?0:den/(den-num);
    assert(num==den || den%(den-num)==0);
}
Rat evaluate_q(const Poly& a,int q) {
    W n=0; int d=a.size()-1;
    for(int i=0;i<=d;++i) n=n*q+a[i];
    return Rat(n,power(q,d));
}
void artin_tate(Row& r) {
    Poly a=r.P2,z;
    for(int i=0;i<r.rho;++i) { bool ok=divide(a,{1,-r.q},z); assert(ok); a=z; }
    auto e=evaluate_q(a,r.q); r.brdisc=Rat(W(r.q)*e.n,e.d);
    assert(r.brdisc.n>0 && r.brdisc.d==1);
}
void count_points(Row& r) {
    std::array<I,6> s{};
    for(int k=1;k<=5;++k) {
        W n=W(k)*r.P2[k];
        for(int j=1;j<k;++j) n+=W(r.P2[j])*s[k-j];
        s[k]=narrow(-n); r.points[k-1]=1+power(r.q,2*k)+s[k];
    }
}
void bsd(Row& r) {
    auto e=evaluate_q(r.quotient,r.q); I cp=1;
    for(auto f:r.fibers) cp*=f.c;
    r.sha_reg=Rat(W(r.q)*e.n*r.tor*r.tor,W(e.d)*cp);
}
void split_values(Row& r) {
    r.disc=Rat(r.brdisc.n,W(r.brdisc.d)*r.sha);
    r.reg=Rat(r.sha_reg.n,W(r.sha_reg.d)*r.sha);
    assert(r.disc.d==1);
}
void sign(Row& r) { r.sign=r.P2.back()>0?1:-1; }
void torsion_screen(Row& r) {
    Rat h;
    for(auto f:r.fibers) {
        Rat c;
        switch(f.type) {
        case 2: c=Rat(1,2); break;
        case 3: c=Rat(2,3); break;
        case 4: c=Rat(4,3); break;
        case 5: c=Rat(3,2); break;
        case 7: c=Rat((f.n/2)*(f.n-f.n/2),f.n); break;
        case 8: c=Rat(f.n+4,4); break;
        }
        h=Rat(W(h.n)*c.d+W(c.n)*f.d*h.d,W(h.d)*c.d);
    }
    r.tor_trivial=h.n<4*h.d;
    assert(!r.tor_trivial || r.tor==1);
}
double cpu() { timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t); return t.tv_sec+1e-9*t.tv_nsec; }
volatile I checksum=0;
template<class F> void measure(std::vector<Row>& rows,const char* label,F fn,int reps) {
    int n=0; for(const auto& r:rows) n+=r.timed;
    double start=cpu();
    for(int k=0;k<reps;++k) for(auto& r:rows) if(r.timed) fn(r);
    double seconds=cpu()-start;
    for(const auto& r:rows) checksum=checksum+r.rank+r.rhog+r.height+r.points[0]+r.nsdegree+r.sign+r.disc.n;
    std::cerr<<"TIMING,"<<rows[0].q<<","<<label<<","<<n<<","<<reps<<","<<seconds<<"\n";
}
void vector_json(const Poly& p) { std::cout<<"["; for(size_t i=0;i<p.size();++i) std::cout<<(i?",":"")<<p[i]; std::cout<<"]"; }
int main(int argc,char**argv) {
    assert(argc>=2); int reps=argc>2?std::stoi(argv[2]):100;
    for(int q:{2,3}) for(auto entry:cyclotomics) {
        auto& c=entry.second; I scale=1;
        for(auto& a:c) { a*=scale; scale*=q; }
        if(c[0]<0) for(auto& a:c) a=-a;
        scaled[q].push_back(entry);
    }
    std::ifstream f(argv[1]); std::vector<Row> rows;
    while(true) {
        Row r; int n; if(!(f>>r.q>>r.key>>r.timed>>r.tor>>r.sha>>r.known>>n)) break;
        r.L.resize(n); for(auto& a:r.L) f>>a;
        f>>n; r.fibers.resize(n); for(auto& x:r.fibers) f>>x.d>>x.type>>x.n>>x.c;
        components(r); p2(r); rational_rank(r); geometric(r); height(r);
        artin_tate(r); count_points(r); bsd(r); split_values(r); sign(r); torsion_screen(r);
        rows.push_back(r);
    }
    assert(!rows.empty());
    measure(rows,"components",components,reps);
    measure(rows,"P2",p2,reps);
    measure(rows,"rho",rational_rank,reps);
    measure(rows,"geometric",geometric,reps);
    measure(rows,"height",height,reps);
    measure(rows,"br_disc",artin_tate,reps);
    measure(rows,"points",count_points,reps);
    measure(rows,"bsd",bsd,reps);
    measure(rows,"split",split_values,reps);
    measure(rows,"sign",sign,reps);
    measure(rows,"torsion_screen",torsion_screen,reps);
    for(auto& r:rows) {
        std::cout<<"{\"q\":"<<r.q<<",\"key\":"<<r.key<<",\"P2\":"; vector_json(r.P2);
        std::cout<<",\"cycles\":"; vector_json(Poly(r.cycles.begin(),r.cycles.end()));
        std::cout<<",\"rank\":"<<r.rank<<",\"geom_rank\":"<<r.geom<<",\"rho\":"<<r.rho<<",\"rho_geom\":"<<r.rhog;
        std::cout<<",\"height\":"<<r.height<<",\"br_disc\":"<<r.brdisc.n<<",\"points\":";
        vector_json(Poly(r.points.begin(),r.points.end()));
        std::cout<<",\"zeta_sign\":"<<r.sign<<",\"ns_degree\":"<<r.nsdegree;
        std::cout<<",\"sha_reg\":\""<<r.sha_reg.str()<<"\",\"sha_candidate\":"<<r.sha<<",\"disc_candidate\":\""<<r.disc.str();
        std::cout<<"\",\"reg_candidate\":\""<<r.reg.str()<<"\",\"split_known\":"<<r.known<<",\"tor_trivial\":"<<r.tor_trivial<<"}\n";
    }
}
