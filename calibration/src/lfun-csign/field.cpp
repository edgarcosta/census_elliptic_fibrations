#include "field.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

static const int polys2[12][13]={
 {1,1},{1,1,1},{1,1,0,1},{1,1,0,0,1},{1,0,1,0,0,1},
 {1,1,0,1,1,0,1},{1,1,0,0,0,0,0,1},{1,0,1,1,1,0,0,0,1},
 {1,0,0,0,1,0,0,0,0,1},{1,1,1,1,0,1,1,0,0,0,1},
 {1,0,1,0,0,0,0,0,0,0,0,1},{1,1,0,1,0,1,1,1,0,0,0,0,1}};
static const int polys3[12][13]={
 {1,1},{2,2,1},{1,2,0,1},{2,0,0,2,1},{1,2,0,0,0,1},
 {2,2,1,0,2,0,1},{1,0,2,0,0,0,0,1},{2,2,2,0,1,2,0,0,1},
 {1,1,2,2,0,0,0,0,0,1},{2,1,0,0,2,2,2,0,0,0,1},
 {1,0,2,0,0,0,0,0,0,0,0,1},{2,0,1,0,1,1,1,0,0,0,0,0,1}};

static void fft3(std::vector<std::complex<double>>& a,bool inverse) {
    const int n=a.size();
    const std::complex<double> w(-0.5,(inverse?1:-1)*std::sqrt(3.0)/2);
    for(int stride=1;stride<n;stride*=3)
        for(int block=0;block<n;block+=3*stride)
            for(int j=0;j<stride;++j) {
                auto x=a[block+j],y=a[block+j+stride],z=a[block+j+2*stride];
                a[block+j]=x+y+z;
                a[block+j+stride]=x+w*y+std::conj(w)*z;
                a[block+j+2*stride]=x+std::conj(w)*y+w*z;
            }
    if(inverse) for(auto& v:a) v/=n;
}

static constexpr int64_t modulus=1000000009;
static int64_t modpow(int64_t a,int64_t n) {
    int64_t r=1;
    for(;n;n>>=1,a=a*a%modulus)if(n&1)r=r*a%modulus;
    return r;
}
static void exact3(std::vector<int64_t>& a,bool inverse) {
    const int n=a.size();
    int64_t w=modpow(2,(modulus-1)/3);
    if(w==1 || (w*w+w+1)%modulus)throw std::runtime_error("invalid cube root");
    if(inverse)w=w*w%modulus;
    const int64_t ww=w*w%modulus;
    for(int stride=1;stride<n;stride*=3)
        for(int block=0;block<n;block+=3*stride)
            for(int j=0;j<stride;++j) {
                auto x=a[block+j],y=a[block+j+stride],z=a[block+j+2*stride];
                a[block+j]=(x+y+z)%modulus;
                a[block+j+stride]=(x+w*y+ww*z)%modulus;
                a[block+j+2*stride]=(x+ww*y+w*z)%modulus;
            }
    if(inverse) {
        int64_t ni=modpow(n,modulus-2);
        for(auto& v:a)v=v*ni%modulus;
    }
}

Field::Field(int characteristic,int degree):p(characteristic),d(degree),Q(1),mask((1<<d)-1) {
    for(int i=0;i<d;++i) Q*=p;
    m=Q-1;
    index_mask.resize(1<<d);
    for(int i=1;i<(1<<d);++i) {
        int b=__builtin_ctz(unsigned(i)),v=1;
        for(int j=0;j<b;++j)v*=3;
        index_mask[i]=index_mask[i&(i-1)]+v;
    }
    exps.resize(m); logs.assign(Q,-1);
    const int* poly=p==2 ? polys2[d-1] : polys3[d-1];
    uint32_t mod=0;
    for(int i=0;i<d;++i) if(poly[i]) mod|=uint32_t(1)<<(i+(p==3&&poly[i]==2?16:0));
    uint32_t a=1;
    for(int e=0;e<m;++e) {
        exps[e]=a;
        if(logs[index(a)]!=-1)throw std::runtime_error("nonprimitive field polynomial");
        logs[index(a)]=e;
        int top=((a>>(d-1))&1)+(p==3 ? 2*((a>>(d+15))&1):0);
        a=p==2 ? (a<<1)&mask : ((a<<1)&mask)|(((a>>16<<1)&mask)<<16);
        if(top) a=add(a,p==2 || top==2 ? mod : neg(mod));
    }
    if(a!=1)throw std::runtime_error("field cycle did not close");
    std::vector<int> tb(d);
    for(int i=0;i<d;++i) {
        uint32_t b=1u<<i,t=0;
        for(int j=0;j<d;++j){t=add(t,b); b=pow(b,p);}
        tb[i]=p==2 ? int(t) : int(t&1)+2*int((t>>16)&1);
        if(tb[i]>=p)throw std::runtime_error("field trace is not constant");
    }
    traces.assign(Q,0);
    for(int e=0;e<m;++e) {
        uint32_t b=exps[e]; int t=0;
        for(int i=0;i<d;++i) t+=tb[i]*(((b>>i)&1)+(p==3?2*((b>>(i+16))&1):0));
        traces[index(b)]=t%p;
    }
    ordinary.resize(m+2);
    if(p==3) {
        std::vector<std::complex<double>> h(Q),ch(Q);
        std::vector<int64_t> eh(Q),ec(Q);
        for(int e=0;e<m;++e) {
            uint32_t b=exps[e];
            h[index(neg(add(pow(b,3),pow(b,2))))]+=1;
            ++eh[index(neg(add(pow(b,3),pow(b,2))))];
            ch[index(b)]=e&1?-1:1;
            ec[index(b)]=e&1?modulus-1:1;
            gauss[trace(b)]+=e&1?-1:1;
        }
        h[0]+=1;
        ++eh[0];
        // Each convolution coefficient is a sum of Q values in {-1,0,1}.
        // M > 2Q makes its centered residue the unique integer coefficient.
        if(modulus<=2LL*Q)throw std::runtime_error("trace modulus too small");
        exact3(eh,false);exact3(ec,false);
        for(int i=0;i<Q;++i)eh[i]=eh[i]*ec[i]%modulus;
        exact3(eh,true);
        fft3(h,false); fft3(ch,false);
        for(int i=0;i<Q;++i) h[i]*=ch[i];
        fft3(h,true);
        for(int i=0;i<Q;++i) {
            if(eh[i]>modulus/2)eh[i]-=modulus;
            if(!std::isfinite(h[i].real()) || !std::isfinite(h[i].imag()) || std::abs(h[i].real()-double(eh[i]))>1e-5 || std::abs(h[i].imag())>1e-5)
                throw std::runtime_error("exact F_3 trace table mismatch at degree "+std::to_string(d)+" index "+std::to_string(i));
        }
        for(int e=0;e<m;++e) {
            auto t=h[index(exps[e])];
            if(std::abs(t.real()-std::round(t.real()))>1e-5 || std::abs(t.imag())>1e-5)
                throw std::runtime_error("nonintegral Fourier trace");
            ordinary[e]=int16_t(-eh[index(exps[e])]);
        }
        if(std::llround(-h[0].real())!=1)throw std::runtime_error("nodal trace is not one");
    } else {
        std::vector<int> h(Q);
        for(int e=0;e<m;++e)h[exp(-2LL*e)]=trace(exps[e])?-1:1;
        for(int stride=1;stride<Q;stride*=2)
            for(int base=0;base<Q;base+=2*stride)
                for(int j=0;j<stride;++j) {
                    int x=h[base+j],y=h[base+j+stride];
                    h[base+j]=x+y; h[base+j+stride]=x-y;
                }
        for(int e=0;e<m;++e) {
            uint32_t b=exps[e]; int dual=0;
            for(int j=0;j<d;++j)dual|=trace(mul(b,1u<<j))<<j;
            ordinary[e]=int16_t(-h[dual]);
        }
        if(-h[0]!=1)throw std::runtime_error("nodal trace is not one");
    }
    for(int e=0;e<m;++e) {
        int f=e,orbit=0,smallest=e;
        do {smallest=std::min(smallest,f); f=int(int64_t(f)*p%m); ++orbit;}while(f!=e);
        if(orbit==d && smallest==e) reps.push_back(exps[e]);
    }
    if(d==1)reps.push_back(0);
    for(auto& col:powers)col.resize(reps.size());
    for(size_t i=0;i<reps.size();++i){
        uint32_t v=1;
        for(int j=0;j<=24;++j){powers[j][i]=v; v=mul(v,reps[i]);}
    }
}

int Field::ss3(uint32_t a,uint32_t b) const {
    if(!a)return 0;
    int e=log(neg(a));
    if(e&1)return 0;
    uint32_t r=exp(e/2);
    return -3*chi(r)*gauss[trace(div(b,pow(r,3)))];
}

int Field::ss2(uint32_t a2,uint32_t a3,uint32_t a4,uint32_t a6) const {
    if(!a3)return 0;
    uint32_t inv=exp(-2LL*log(a3));
    uint32_t u=inv,v=mul(a2,inv),w=mul(a4,inv),z=mul(a6,inv);
    unsigned lin=0,active=(1u<<d)-1;
    unsigned cross[12]{};
    int c=trace(z),factor=1;
    for(int i=0;i<d;++i) {
        uint32_t x=1u<<i;
        lin|=unsigned(trace(add(add(mul(u,pow(x,3)),mul(v,pow(x,2))),mul(w,x))))<<i;
        for(int j=0;j<i;++j) {
            uint32_t y=1u<<j;
            int b=trace(mul(u,add(mul(pow(x,2),y),mul(x,pow(y,2)))));
            if(b){cross[i]|=1u<<j;cross[j]|=1u<<i;}
        }
    }
    while(true) {
        int i=-1,j=-1;
        for(int k=0;k<d;++k)if((active>>k&1)&&(cross[k]&active)) {i=k;j=__builtin_ctz(cross[k]&active);break;}
        if(i<0)break;
        active&=~((1u<<i)|(1u<<j));
        unsigned A=cross[i]&active,B=cross[j]&active;
        int ai=(lin>>i)&1,bj=(lin>>j)&1;
        c^=ai&bj;
        lin^=(ai?B:0)^(bj?A:0)^(A&B);
        for(int k=0;k<d;++k)if(active>>k&1) {
            unsigned row=((A>>k&1)?B:0)^((B>>k&1)?A:0);
            cross[k]^=row&~(1u<<k);
        }
        factor*=2;
    }
    if(lin&active)return 0;
    return -(c?-1:1)*factor*(1<<__builtin_popcount(active));
}
