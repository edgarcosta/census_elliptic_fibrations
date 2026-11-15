#ifndef LFUN_FIELD_HPP
#define LFUN_FIELD_HPP
#include <array>
#include <complex>
#include <cstdint>
#include <vector>

struct Field {
    int p,d,Q,m,mask;
    std::vector<uint32_t> exps, index_mask, reps;
    std::vector<int32_t> logs;
    std::vector<int16_t> ordinary;
    std::vector<uint8_t> traces;
    std::array<std::vector<uint32_t>,25> powers;
    int gauss[3]{};
    Field(int characteristic,int degree);
    uint32_t add(uint32_t a,uint32_t b) const {
        if(p==2) return a^b;
        uint32_t al=a&65535,ah=a>>16,bl=b&65535,bh=b>>16;
        uint32_t an=al|ah,bn=bl|bh;
        return ((al&~bn)|(bl&~an)|(ah&bh)) | (((ah&~bn)|(bh&~an)|(al&bl))<<16);
    }
    uint32_t neg(uint32_t a) const {return p==2 ? a : (a>>16)|(a<<16);}
    uint32_t sub(uint32_t a,uint32_t b) const {return add(a,neg(b));}
    int index(uint32_t a) const {return p==2 ? int(a) : int(index_mask[a&mask]+2*index_mask[a>>16]);}
    int log(uint32_t a) const {return logs[index(a)];}
    uint32_t exp(int64_t e) const {e%=m; if(e<0)e+=m; return exps[e];}
    uint32_t mul(uint32_t a,uint32_t b) const {return a&&b ? exp(log(a)+log(b)) : 0;}
    uint32_t div(uint32_t a,uint32_t b) const {return a ? exp(log(a)-log(b)) : 0;}
    uint32_t pow(uint32_t a,int64_t e) const {return a ? exp(log(a)*e) : 0;}
    int trace(uint32_t a) const {return traces[index(a)];}
    int chi(uint32_t a) const {return a ? (log(a)&1 ? -1 : 1) : 0;}
    int ss3(uint32_t a,uint32_t b) const;
    int ss2(uint32_t a2,uint32_t a3,uint32_t a4,uint32_t a6) const;
};
#endif
