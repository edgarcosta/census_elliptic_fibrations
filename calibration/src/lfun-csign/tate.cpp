#include "tate.hpp"
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEX.h>
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
using namespace NTL;
namespace {
using P=zz_pX;
using E=zz_pE;
using Model=std::array<P,5>;
P poly(const Poly& a) {
    P f;
    for(int i=0;i<=24;++i)SetCoeff(f,i,a[i]);
    return f;
}
P quotient(const P& a,const P& pi,int n=1) {
    P b=a,q,r;
    while(n--) {
        DivRem(q,r,b,pi);
        if(!IsZero(r))throw std::runtime_error("inexact uniformizer division");
        b=q;
    }
    return b;
}
int valuation(P a,const P& pi) {
    if(IsZero(a))return std::numeric_limits<int>::max();
    int v=0;
    while(IsZero(a%pi)){a=quotient(a,pi);++v;}
    return v;
}
void transform(Model& a,const P& r,const P& s,const P& t) {
    auto a1=a[0],a2=a[1],a3=a[2],a4=a[3],a6=a[4];
    a[4]=a6+r*(a4+r*(a2+r))-t*(a3+r*a1+t);
    a[3]=a4-s*a3+2*r*a2-(t+r*s)*a1+3*r*r-2*s*t;
    a[2]=a3+r*a1+2*t;
    a[1]=a2-s*a1+3*r-s*s;
    a[0]=a1+2*s;
}
struct Invariants {
    P b2,b4,b6,b8;
    explicit Invariants(const Model& a) {
        const auto &a1=a[0],&a2=a[1],&a3=a[2],&a4=a[3],&a6=a[4];
        b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;
        b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4;
    }
    P delta() const {return -b2*b2*b8-8*b4*b4*b4-27*b6*b6+9*b2*b4*b6;}
};
struct Residue {
    int p,d;
    long Q;
    explicit Residue(const P& pi):p(zz_p::modulus()),d(deg(pi)),Q(1) {
        zz_pE::init(pi);
        for(int i=0;i<d;++i)Q*=p;
    }
    E reduce(const P& a) const {return conv<E>(a);}
    P lift(const E& a) const {return rep(a);}
    P mod(const P& a) const {return lift(reduce(a));}
    P inverse(const P& a) const {return lift(inv(reduce(a)));}
    P root(const P& a) const {return lift(power(reduce(a),Q/p));}
    int character(const E& a) const {
        E v=power(a,(Q-1)/2);
        return IsOne(v)?1:(IsZero(v)?0:-1);
    }
    int trace(E a) const {
        E sum;
        for(int i=0;i<d;++i){sum+=a;a=power(a,p);}
        return rep(coeff(rep(sum),0));
    }
    bool quadratic(E a,E b,E c) const {
        if(IsZero(a))return !IsZero(b)||IsZero(c);
        if(p==3)return character(b*b-4*a*c)>=0;
        if(IsZero(b))return true;
        return trace(a*c/(b*b))==0;
    }
    bool quadratic(const P& a,const P& b,const P& c) const {
        return quadratic(reduce(a),reduce(b),reduce(c));
    }
    int cubicroots(const P& b,const P& c,const P& d0) const {
        zz_pEX f,x;
        SetCoeff(f,3,1);SetCoeff(f,2,reduce(b));SetCoeff(f,1,reduce(c));SetCoeff(f,0,reduce(d0));
        SetCoeff(x,1,1);
        return deg(GCD(f,PowerXMod(Q,zz_pEXModulus(f))-x));
    }
};
LocalData local(Model a,const P& pi,int vd,const std::string& place) {
    Residue k(pi);
    const int p=k.p;
    int scalings=0;
    auto finish=[&](std::string symbol,int m,int cp,int trace=0) {
        int conductor=vd+1-m;
        if(conductor<0)throw std::runtime_error("negative Ogg conductor");
        return LocalData{place,symbol,int(deg(pi)),vd,conductor,cp,trace,scalings};
    };
    for(;;) {
        if(vd==0)return finish("I0",1,1);
        Invariants b(a);
        E primary=k.reduce(p==2?a[0]:b.b2);
        if(!IsZero(primary)) {
            int tr;
            if(p==3)tr=k.character(primary);
            else {
                E a1=primary;
                tr=k.trace(k.reduce(a[1])/(a1*a1)+k.reduce(a[2])/(a1*a1*a1))?-1:1;
            }
            return finish("I"+std::to_string(vd),vd,tr==1?vd:(vd%2?1:2),tr);
        }
        P r,t;
        if(p==2) {
            r=k.root(a[3]);
            t=k.root(((r+a[1])*r+a[3])*r+a[4]);
        } else {
            r=k.root(-b.b6);
            t=a[0]*r+a[2];
        }
        r=k.mod(r);t=k.mod(t);
        transform(a,r,P(),t);
        b=Invariants(a);
        if(valuation(a[4],pi)<2)return finish("II",1,1);
        if(valuation(b.b8,pi)<3)return finish("III",2,2);
        if(valuation(b.b6,pi)<3)
            return finish("IV",3,k.quadratic(P(1),quotient(a[2],pi),-quotient(a[4],pi,2))?3:1);
        P s;
        if(p==2) {
            s=k.root(a[1]);
            t=pi*k.root(quotient(a[4],pi,2));
        } else {s=a[0];t=a[2];}
        transform(a,P(),s,t);
        P cb=quotient(a[1],pi),cc=quotient(a[3],pi,2),cd=quotient(a[4],pi,3);
        E rb=k.reduce(cb),rc=k.reduce(cc),rd=k.reduce(cd);
        E w=27*rd*rd-rb*rb*rc*rc+4*rb*rb*rb*rd-18*rb*rc*rd+4*rc*rc*rc;
        E x=3*rc-rb*rb;
        if(!IsZero(w))return finish("I0*",5,1+k.cubicroots(cb,cc,cd));
        if(!IsZero(x)) {
            r=pi*(p==2?k.root(cc):k.mod(cc*k.inverse(cb)));
            transform(a,r,P(),P());
            int ix=3,iy=3;
            P mx=pi*pi,my=mx;
            for(;;) {
                P a2t=quotient(a[1],pi),a3t=quotient(a[2],my),a4t=quotient(a[3],pi*mx),a6t=quotient(a[4],mx*my);
                E u=k.reduce(a3t),v=k.reduce(a6t);
                int cp;
                if(IsZero(u*u+4*v)) {
                    t=my*(p==2?k.root(a6t):k.mod(a3t));
                    transform(a,P(),P(),t);
                    my*=pi;++iy;
                    a2t=quotient(a[1],pi);a4t=quotient(a[3],pi*mx);a6t=quotient(a[4],mx*my);
                    E aa=k.reduce(a2t),bb=k.reduce(a4t),ccc=k.reduce(a6t);
                    if(IsZero(bb*bb-4*ccc*aa)) {
                        r=mx*(p==2?k.root(a6t*k.inverse(a2t)):k.mod(a4t*k.inverse(a2t)));
                        transform(a,r,P(),P());
                        mx*=pi;++ix;
                        if(ix+iy>vd+5)throw std::runtime_error("Tate I* loop exceeded discriminant");
                        continue;
                    }
                    cp=k.quadratic(aa,bb,ccc)?4:2;
                } else cp=k.quadratic(E(1),u,-v)?4:2;
                int n=ix+iy-5;
                return finish("I"+std::to_string(n)+"*",n+5,cp);
            }
        }
        r=pi*(p==2?k.mod(cb):k.root(-cd));
        transform(a,r,P(),P());
        P a3t=quotient(a[2],pi,2),a6t=quotient(a[4],pi,4);
        E u=k.reduce(a3t),v=k.reduce(a6t);
        if(!IsZero(u*u+4*v))return finish("IV*",7,k.quadratic(E(1),u,-v)?3:1);
        t=pi*pi*(p==2?-k.root(a6t):k.mod(a3t));
        transform(a,P(),P(),t);
        if(valuation(a[3],pi)<4)return finish("III*",8,2);
        if(valuation(a[4],pi)<6)return finish("II*",9,1);
        constexpr int weight[]={1,2,3,4,6};
        for(int j=0;j<5;++j)a[j]=quotient(a[j],pi,weight[j]);
        vd-=12;++scalings;
        if(vd<0)throw std::runtime_error("Tate scaling exceeded discriminant");
    }
}
uint64_t encoding(const P& f,int p) {
    uint64_t n=0;
    for(int i=deg(f);i>=0;--i)n=n*p+rep(coeff(f,i));
    return n;
}
}
std::vector<LocalData> tate(const Surface& s,bool* everywhere_minimal,SignData* sign) {
    if(everywhere_minimal)*everywhere_minimal=true;
    zz_p::init(s.p);
    Model a;
    for(int j=0;j<5;++j)a[j]=poly(s.a[j]);
    auto frobenius=[&](const Model& original) {
        Model twist;
        for(int j=0;j<5;++j)for(int i=0;i<=deg(original[j]);++i)
            SetCoeff(twist[j],s.p*i,coeff(original[j],i));
        return twist;
    };
    Model twist;
    if(sign){*sign=SignData{};twist=frobenius(a);}
    int sum0=0,sum1=0;
    auto compare=[&](const LocalData& original,const LocalData& twisted) {
        auto vp=[&](int n){int v=0;while(n%s.p==0){n/=s.p;++v;}return v;};
        sum0+=original.degree*original.vd;
        sum1+=twisted.degree*twisted.vd;
        sign->tamagawa_exponent+=vp(twisted.tamagawa)-vp(original.tamagawa);
    };
    P delta=Invariants(a).delta();
    if(IsZero(delta))throw std::runtime_error("singular generic fibre for key "+std::to_string(s.key));
    int infinity_vd=24-deg(delta);
    MakeMonic(delta);
    vec_pair_zz_pX_long factors;
    CanZass(factors,delta);
    std::vector<LocalData> result;
    for(long i=0;i<factors.length();++i) {
        const P& pi=factors[i].a;
        auto data=local(a,pi,factors[i].b,std::to_string(encoding(pi,s.p)));
        if(sign)compare(data,local(twist,pi,s.p*factors[i].b,data.place));
        if(everywhere_minimal&&data.scalings)*everywhere_minimal=false;
        if(data.conductor)result.push_back(data);
    }
    if(infinity_vd>0) {
        Model chart;
        constexpr int weight[]={2,4,6,8,12};
        for(int j=0;j<5;++j)for(int i=0;i<=weight[j];++i)SetCoeff(chart[j],weight[j]-i,s.a[j][i]);
        P pi;SetCoeff(pi,1,1);
        auto data=local(chart,pi,infinity_vd,"inf");
        // The twist chart has weights 2*p*i, before local minimization.
        if(sign)compare(data,local(frobenius(chart),pi,s.p*infinity_vd,"inf"));
        if(everywhere_minimal&&data.scalings)*everywhere_minimal=false;
        if(data.conductor)result.push_back(data);
    }
    if(sign) {
        if(sum0%12||sum1%12)throw std::runtime_error("nonintegral elliptic surface height");
        sign->height=sum0/12;sign->twist_height=sum1/12;
        int exponent=sign->twist_height-sign->height+sign->tamagawa_exponent;
        sign->sign=exponent%2?-1:1;
    }
    std::sort(result.begin(),result.end(),[](const LocalData& a0,const LocalData& b0){return a0.place<b0.place;});
    return result;
}
