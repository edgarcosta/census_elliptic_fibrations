#include "counter.hpp"
#ifdef USE_AVX2
#include <immintrin.h>
static inline __m256i add3(__m256i a,__m256i b) {
    auto mask=_mm256_set1_epi32(65535);
    auto al=_mm256_and_si256(a,mask),ah=_mm256_srli_epi32(a,16);
    auto bl=_mm256_and_si256(b,mask),bh=_mm256_srli_epi32(b,16);
    auto an=_mm256_or_si256(al,ah),bn=_mm256_or_si256(bl,bh);
    auto lo=_mm256_or_si256(_mm256_or_si256(_mm256_andnot_si256(bn,al),_mm256_andnot_si256(an,bl)),_mm256_and_si256(ah,bh));
    auto hi=_mm256_or_si256(_mm256_or_si256(_mm256_andnot_si256(bn,ah),_mm256_andnot_si256(an,bh)),_mm256_and_si256(al,bl));
    return _mm256_or_si256(lo,_mm256_slli_epi32(hi,16));
}
#endif

void evaluate(const Field& f,const Poly& poly,std::vector<uint32_t>& out) {
    const int n=f.reps.size();
    out.assign(n,0);
    for(int j=0;j<=24;++j)if(poly[j]) {
        const auto* in=f.powers[j].data();
        int i=0;
#ifdef USE_AVX2
        for(;i+8<=n;i+=8) {
            auto v=_mm256_loadu_si256((const __m256i*)(in+i));
            if(poly[j]==2)v=_mm256_or_si256(_mm256_srli_epi32(v,16),_mm256_slli_epi32(v,16));
            auto a=_mm256_loadu_si256((const __m256i*)(out.data()+i));
            a=f.p==2 ? _mm256_xor_si256(a,v) : add3(a,v);
            _mm256_storeu_si256((__m256i*)(out.data()+i),a);
        }
#endif
        for(;i<n;++i)out[i]=f.add(out[i],poly[j]==2?f.neg(in[i]):in[i]);
    }
}

int64_t sum_ordinary3(const Field& f,const std::vector<uint32_t>& c,const std::vector<uint32_t>& delta,bool& all_ordinary) {
    const int n=c.size(),m=f.m;
    all_ordinary=true;
    int64_t sum=0; int i=0;
#ifdef USE_AVX2
    auto mm=_mm256_set1_epi32(m),zero=_mm256_setzero_si256();
    auto one=_mm256_set1_epi32(1),six=_mm256_set1_epi32(6),sixm=_mm256_set1_epi32(6*m);
    auto mask=_mm256_set1_epi32(f.mask);
    auto recip=_mm256_set1_ps(1.0f/m);
    for(;i+8<=n;i+=8) {
        auto a=_mm256_loadu_si256((const __m256i*)(c.data()+i));
        if(_mm256_movemask_epi8(_mm256_cmpeq_epi32(a,zero))) {all_ordinary=false;return 0;}
        auto b=_mm256_loadu_si256((const __m256i*)(delta.data()+i));
        auto idxa=_mm256_add_epi32(_mm256_i32gather_epi32((const int*)f.index_mask.data(),_mm256_and_si256(a,mask),4),
            _mm256_slli_epi32(_mm256_i32gather_epi32((const int*)f.index_mask.data(),_mm256_srli_epi32(a,16),4),1));
        auto idxb=_mm256_add_epi32(_mm256_i32gather_epi32((const int*)f.index_mask.data(),_mm256_and_si256(b,mask),4),
            _mm256_slli_epi32(_mm256_i32gather_epi32((const int*)f.index_mask.data(),_mm256_srli_epi32(b,16),4),1));
        auto la=_mm256_i32gather_epi32(f.logs.data(),idxa,4),lb=_mm256_i32gather_epi32(f.logs.data(),idxb,4);
        auto e=_mm256_add_epi32(_mm256_sub_epi32(lb,_mm256_mullo_epi32(la,six)),sixm);
        auto quotient=_mm256_cvttps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(e),recip));
        e=_mm256_sub_epi32(e,_mm256_mullo_epi32(quotient,mm));
        e=_mm256_add_epi32(e,_mm256_and_si256(_mm256_cmpgt_epi32(zero,e),mm));
        e=_mm256_sub_epi32(e,_mm256_and_si256(_mm256_cmpgt_epi32(e,_mm256_sub_epi32(mm,one)),mm));
        auto tr=_mm256_i32gather_epi32((const int*)f.ordinary.data(),e,2);
        tr=_mm256_srai_epi32(_mm256_slli_epi32(tr,16),16);
        tr=_mm256_blendv_epi8(tr,one,_mm256_cmpeq_epi32(b,zero));
        auto sign=_mm256_sub_epi32(zero,_mm256_and_si256(la,one));
        tr=_mm256_sub_epi32(_mm256_xor_si256(tr,sign),sign);
        alignas(32)int vals[8]; _mm256_store_si256((__m256i*)vals,tr);
        for(int v:vals)sum+=v;
    }
#endif
    for(;i<n;++i) {
        if(!c[i]){all_ordinary=false;return 0;}
        int a=f.chi(c[i]);
        if(delta[i])a*=f.ordinary[(f.log(delta[i])-6LL*f.log(c[i])%m+m)%m];
        sum+=a;
    }
    return sum;
}
