#include "counter.hpp"
#include "tate.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

static Poly plus(const Poly& a,const Poly& b,int p,int sign=1) {
    Poly c{};
    for(int i=0;i<=24;++i)c[i]=(a[i]+sign*b[i]+p)%p;
    return c;
}
static Poly times(const Poly& a,const Poly& b,int p) {
    Poly c{};
    for(int i=0;i<=24;++i)if(a[i])for(int j=0;i+j<=24;++j)c[i+j]=(c[i+j]+a[i]*b[j])%p;
    return c;
}
void Surface::invariants() {
    auto square=[&](const Poly& v){return times(v,v,p);};
    if(p==3) {
        primary=plus(a[1],square(a[0]),p);
        A=plus(a[3],times(a[0],a[2],p),p,-1);
        B=plus(a[4],square(a[2]),p);
        delta=plus(plus(times(square(primary),times(primary,B,p),p),times(square(A),A,p),p),times(square(primary),square(A),p),p,-1);
    } else {
        primary=a[0];
        auto b2=square(a[0]),b4=times(a[0],a[2],p),b6=square(a[2]);
        auto b8=plus(plus(times(b2,a[4],p),times(b4,a[3],p),p),plus(times(a[1],b6,p),square(a[3]),p),p);
        delta=plus(plus(times(square(b2),b8,p),square(b6),p),times(times(b2,b4,p),b6,p),p);
    }
}
static uint32_t value(const Field& f,const Poly& a,int i) {
    uint32_t v=0;
    for(int j=0;j<=24;++j)if(a[j])v=f.add(v,a[j]==2?f.neg(f.powers[j][i]):f.powers[j][i]);
    return v;
}
static int trace_at(const Surface& s,const Field& f,int i,uint32_t c,uint32_t delta) {
    if(s.p==3) {
        if(c) {
            if(!delta)return f.chi(c);
            return f.chi(c)*f.ordinary[(f.log(delta)-6LL*f.log(c)%f.m+f.m)%f.m];
        }
        return f.ss3(value(f,s.A,i),value(f,s.B,i));
    }
    uint32_t a2=value(f,s.a[1],i),a3=value(f,s.a[2],i);
    if(c) {
        uint32_t twist=f.add(f.div(a2,f.pow(c,2)),f.div(a3,f.pow(c,3)));
        int sign=f.trace(twist)?-1:1;
        if(!delta)return sign;
        return sign*f.ordinary[(f.log(delta)-12LL*f.log(c)%f.m+f.m)%f.m];
    }
    return f.ss2(a2,a3,value(f,s.a[3],i),value(f,s.a[4],i));
}
static void local_series(int a,bool good,int qd,int d,int limit,int64_t* S) {
    int64_t prev=2,cur=a;
    for(int r=1;r*d<=limit;++r) {
        S[r*d]+=d*cur;
        int64_t next=good ? a*cur-qd*prev : a*cur;
        prev=cur;cur=next;
    }
}
static void infinity(const Surface& s,const Field& f,int limit,int64_t* S) {
    constexpr int weight[]={2,4,6,8,12};
    uint32_t a[5];
    for(int j=0;j<5;++j){int c=s.a[j][weight[j]];a[j]=c==2?f.neg(1):c;}
    int tr; uint32_t delta;
    if(s.p==3) {
        uint32_t c=f.add(a[1],f.pow(a[0],2));
        uint32_t A=f.sub(a[3],f.mul(a[0],a[2])),B=f.add(a[4],f.pow(a[2],2));
        delta=f.sub(f.add(f.pow(A,3),f.mul(f.pow(c,3),B)),f.mul(f.pow(c,2),f.pow(A,2)));
        tr=c ? (delta ? f.chi(c)*f.ordinary[(f.log(delta)-6LL*f.log(c)%f.m+f.m)%f.m] : f.chi(c)) : f.ss3(A,B);
    } else {
        uint32_t b2=f.pow(a[0],2),b4=f.mul(a[0],a[2]),b6=f.pow(a[2],2);
        uint32_t b8=f.add(f.add(f.mul(b2,a[4]),f.mul(b4,a[3])),f.add(f.mul(a[1],b6),f.pow(a[3],2)));
        delta=f.add(f.add(f.mul(f.pow(b2,2),b8),f.pow(b6,2)),f.mul(f.mul(b2,b4),b6));
        if(a[0]) {
            int sign=f.trace(f.add(f.div(a[1],f.pow(a[0],2)),f.div(a[2],f.pow(a[0],3))))?-1:1;
            tr=delta ? sign*f.ordinary[(f.log(delta)-12LL*f.log(a[0])%f.m+f.m)%f.m] : sign;
        } else tr=f.ss2(a[1],a[2],a[3],a[4]);
    }
    local_series(tr,delta!=0,s.p,1,limit,S);
}

using Clock=std::chrono::steady_clock;
static double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static int64_t ipow(int a,int n){int64_t v=1;while(n--)v*=a;return v;}
static std::string coeffs(const std::vector<int64_t>& c){std::ostringstream out;for(size_t j=0;j<c.size();++j){if(j)out<<' ';out<<c[j];}return out.str();}
static std::vector<int64_t> complete(const std::vector<int64_t>& c,int p,int D,int sign) {
    std::vector<int64_t> a(D+1);
    for(int i=0;i<=D/2;++i){a[i]=c[i];a[D-i]=sign*ipow(p,D-2*i)*c[i];}
    return a;
}

int main(int argc,char** argv) try {
    if(argc==2 && std::string(argv[1])=="--certify-tables") {
        for(int p:{2,3})for(int d=1;d<=12;++d) {
            Field f(p,d);
            std::cout<<"TABLE_OK|"<<p<<"|"<<d<<"|"<<f.Q<<"|"
                     <<(p==3?"exact radix-3 convolution M=1000000009":"integer Walsh transform")<<'\n';
        }
        return 0;
    }
    int maxdegree=10;
    std::string filename,local_filename;
    bool local_only=false,sign_only=false;
    for(int i=1;i<argc;++i) {
        if(std::string(argv[i])=="--max-degree"&&i+1<argc)maxdegree=std::stoi(argv[++i]);
        else if(std::string(argv[i])=="--local-only")local_only=true;
        else if(std::string(argv[i])=="--sign-only")sign_only=true;
        else if(std::string(argv[i])=="--local-output"&&i+1<argc)local_filename=argv[++i];
        else if(std::string(argv[i])=="--help") {
            std::cout<<"Usage: lfun-avx2 [--max-degree 10|11|12] [--local-only|--sign-only] [--local-output FILE] [input.txt]\nInput: q key sign a1 a2 a3 a4 a6 (base-q polynomial encodings).\n";return 0;
        } else filename=argv[i];
    }
    if(maxdegree<10||maxdegree>12)throw std::runtime_error("max-degree must be 10, 11 or 12");
    if(local_only&&sign_only)throw std::runtime_error("choose local-only or sign-only");
    std::ifstream file;
    if(!filename.empty()){file.open(filename);if(!file)throw std::runtime_error("cannot open input");}
    std::istream& in=filename.empty()?std::cin:file;
    std::array<std::array<std::unique_ptr<Field>,13>,4> fields;
    std::vector<Surface> batch;
    std::string line;
    while(std::getline(in,line)) {
        if(line.empty()||line[0]=='#')continue;
        std::istringstream row(line);
        int p,sign;uint64_t key,codes[5];
        if(!(row>>p>>key>>sign>>codes[0]>>codes[1]>>codes[2]>>codes[3]>>codes[4]))throw std::runtime_error("invalid input row");
        std::string extra;
        if(row>>extra)throw std::runtime_error("extra input fields");
        if((p!=2&&p!=3)||sign< -1||sign>1)throw std::runtime_error("invalid q or sign");
        Surface s(p);s.key=key;s.D=0;s.sign=sign;
        constexpr int weight[]={2,4,6,8,12};
        for(int j=0;j<5;++j){
            s.a[j]={};
            for(int k=0;k<=weight[j];++k){s.a[j][k]=codes[j]%p;codes[j]/=p;}
            if(codes[j])throw std::runtime_error("a-invariant exceeds K3 degree bound");
        }
        batch.push_back(s);
    }
    if(!local_only&&!sign_only)for(int p:{2,3}) {
        if(std::none_of(batch.begin(),batch.end(),[p](const Surface& s){return s.p==p;}))continue;
        auto start=Clock::now();
        for(int d=1;d<=maxdegree;++d)fields[p][d]=std::make_unique<Field>(p,d);
        std::cerr<<"precompute q="<<p<<" degrees="<<maxdegree<<" wall_seconds="<<seconds(start)<<" places10=";
        size_t places=1;for(int d=1;d<=10;++d)places+=fields[p][d]->reps.size();
        size_t bytes=0;
        for(int d=1;d<=maxdegree;++d) {
            const auto& f=*fields[p][d];
            bytes+=4*(f.exps.size()+f.index_mask.size()+f.reps.size()+f.logs.size())+2*f.ordinary.size()+f.traces.size();
            for(const auto& col:f.powers)bytes+=4*col.size();
        }
        std::cerr<<places<<" table_bytes="<<bytes<<'\n';
    }
    std::ofstream local_file;
    if(!local_filename.empty()) {
        local_file.open(local_filename);
        if(!local_file)throw std::runtime_error("cannot open local output");
    }
    std::ostream* local_out=local_only?&std::cout:(!local_filename.empty()?&local_file:nullptr);
    if(local_out)*local_out<<"q,key,place,degree,vd,conductor,tamagawa,kodaira,trace\n";
    if(sign_only)std::cout<<"q,key,sign,height,twist_height,tamagawa_exponent,seconds\n";
    else if(!local_only)std::cout<<"q,key,degree,sign,counted_degree,status,seconds,L_polynomial,other_candidate\n";
    for(auto& s:batch) {
        auto start=Clock::now();
        bool minimal;
        SignData sign_data;
        auto locals=tate(s,&minimal,local_only?nullptr:&sign_data);
        if(!local_only&&!minimal)throw std::runtime_error("counter requires an everywhere minimal K3 model");
        if(!local_only&&sign_data.height!=2)throw std::runtime_error("counter requires K3 height 2");
        if(!local_only&&s.sign&&s.sign!=sign_data.sign)
            throw std::runtime_error("supplied sign conflicts with equation for key "+std::to_string(s.key));
        s.D=-4;
        for(const auto& v:locals) {
            s.D+=v.degree*v.conductor;
        }
        auto write_locals=[&]() {
            if(local_out)for(const auto& v:locals)
                *local_out<<s.p<<','<<s.key<<','<<v.place<<','<<v.degree<<','<<v.vd<<','
                          <<v.conductor<<','<<v.tamagawa<<','<<v.kodaira<<','<<v.trace<<'\n';
        };
        if(local_only){write_locals();continue;}
        if(sign_only) {
            double elapsed=seconds(start);
            write_locals();
            std::cout<<s.p<<','<<s.key<<','<<sign_data.sign<<','<<sign_data.height<<','
                     <<sign_data.twist_height<<','<<sign_data.tamagawa_exponent<<','<<elapsed<<'\n';
            continue;
        }
        if(s.D<0||s.D>20)throw std::runtime_error("conductor outside nonconstant K3 range");
        s.invariants();
        const bool supersingular3=s.p==3&&std::all_of(s.primary.begin(),s.primary.end(),[](int v){return v==0;});
        int64_t S[13]{};
        std::vector<int64_t> c(13);c[0]=1;
        const int sign=sign_data.sign;
        int last=0;
        std::vector<uint32_t> primary,delta;
        infinity(s,*fields[s.p][1],maxdegree,S);
        for(int d=1;d<=std::min(maxdegree,s.D);++d) {
            const auto& f=*fields[s.p][d];
            if(supersingular3) {
                evaluate(f,s.A,primary);evaluate(f,s.B,delta);
                for(size_t i=0;i<f.reps.size();++i)
                    local_series(f.ss3(primary[i],delta[i]),primary[i]!=0,f.Q,d,maxdegree,S);
            } else {
                evaluate(f,s.primary,primary);evaluate(f,s.delta,delta);
                bool all=false;
                if(s.p==3&&d>maxdegree/2)S[d]+=d*sum_ordinary3(f,primary,delta,all);
                if(!all)for(size_t i=0;i<f.reps.size();++i) {
                    int a=trace_at(s,f,i,primary[i],delta[i]);
                    local_series(a,delta[i]!=0,f.Q,d,maxdegree,S);
                }
            }
            __int128 sum=0;
            for(int j=1;j<=d;++j)sum+=(__int128)S[j]*c[d-j];
            if(sum%d)throw std::runtime_error("nonintegral Euler coefficient for key "+std::to_string(s.key));
            c[d]=int64_t(sum/d);last=d;
            if(d*2>=s.D) {
                int i=s.D-d;
                int64_t mirror=ipow(s.p,2*d-s.D)*c[i];
                if(mirror) {
                    int decided=c[d]==mirror?1:(c[d]==-mirror?-1:0);
                    if(!decided)throw std::runtime_error("functional equation failed for key "+std::to_string(s.key));
                    if(sign!=decided)throw std::runtime_error("equation sign conflicts with traces for key "+std::to_string(s.key));
                } else if(c[d])throw std::runtime_error("functional equation zero failed");
                break;
            }
        }
        auto L=complete(c,s.p,s.D,sign);
        double elapsed=seconds(start);
        write_locals();
        std::cout<<s.p<<','<<s.key<<','<<s.D<<','<<sign<<','<<last<<",ok,";
        char time[40];std::snprintf(time,sizeof(time),"%.9f",elapsed);
        std::cout<<time<<','<<coeffs(L)<<",\n";
    }
    return 0;
} catch(const std::exception& e){std::cerr<<"lfun: "<<e.what()<<'\n';return 1;}
