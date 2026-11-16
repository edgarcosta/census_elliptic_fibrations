#include "model.hpp"

#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEX.h>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace census {
namespace {
using NTL::zz_pX;
using Coefficients = std::array<zz_pX, 5>;

struct WeierstrassInvariants {
    zz_pX b2, b4, b6, b8;

    explicit WeierstrassInvariants(const Coefficients& a) {
        const auto& [a1, a2, a3, a4, a6] = a;
        b2 = a1*a1 + 4*a2;
        b4 = a1*a3 + 2*a4;
        b6 = a3*a3 + 4*a6;
        b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4;
    }

    zz_pX discriminant() const {
        return -b2*b2*b8 - 8*b4*b4*b4 - 27*b6*b6 + 9*b2*b4*b6;
    }
};

zz_pX polynomial(const std::vector<int>& coefficients) {
    zz_pX result;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        NTL::SetCoeff(result, i, coefficients[i]);
    }
    return result;
}

bool form_power(const zz_pX& f, int exponent) {
    for (long i = 0; i <= NTL::deg(f); ++i) {
        if (i % exponent != 0 && !NTL::IsZero(NTL::coeff(f, i))) return false;
    }
    return true;
}

std::string fibre_kind(const Coefficients& a, unsigned q) {
    const auto& [a1, a2, a3, a4, a6] = a;
    const zz_pX discriminant = WeierstrassInvariants(a).discriminant();
    bool cusp;
    bool regular;
    if (q == 2) {
        cusp = NTL::IsZero(a1) && NTL::IsZero(a3);
        regular = !form_power(a4, 2) || !form_power(a6, 2);
    } else {
        cusp = NTL::IsZero(a2) && NTL::IsZero(a4);
        regular = !form_power(a6, 3);
    }
    if (!NTL::IsZero(discriminant)) return "smooth";
    return cusp && regular ? "regular cusp" : "nonregular";
}

std::vector<zz_pX> residues(unsigned q, int width) {
    unsigned count = 1;
    for (int i = 0; i < width; ++i) count *= q;
    std::vector<zz_pX> result;
    result.reserve(count);
    for (unsigned code = 0; code < count; ++code) {
        zz_pX f;
        unsigned digits = code;
        for (int i = 0; i < width; ++i) {
            NTL::SetCoeff(f, i, digits % q);
            digits /= q;
        }
        result.push_back(std::move(f));
    }
    return result;
}

bool divisible(const zz_pX& f, const zz_pX& modulus) {
    return NTL::IsZero(f % modulus);
}

bool nonminimal_at(const Coefficients& a, const zz_pX& pi, unsigned q) {
    const auto& [a1, a2, a3, a4, a6] = a;
    const int degree = NTL::deg(pi);
    const zz_pX pi2 = NTL::sqr(pi);
    const zz_pX pi4 = NTL::sqr(pi2);
    const zz_pX pi6 = pi4*pi2;
    if (q == 2) {
        if (!divisible(a1, pi)) return false;
        const zz_pX pi3 = pi2*pi;
        const auto shifts = residues(q, 3*degree);
        for (const auto& slope : residues(q, degree)) {
            const zz_pX r = slope*a1 + NTL::sqr(slope);
            if (!divisible(a3 + r*a1, pi3)) continue;
            for (const auto& v : shifts) {
                const zz_pX c4 = slope*a3 + (v + r*slope)*a1 + NTL::sqr(r);
                if (!divisible(a4 + c4, pi4)) continue;
                const zz_pX c6 = NTL::power(r, 3) + v*a3 + NTL::sqr(v) + r*v*a1;
                if (divisible(a6 + r*a4 + c6, pi6)) return true;
            }
        }
    } else {
        if (!divisible(a2, pi2)) return false;
        for (const auto& r : residues(q, 2*degree)) {
            if (!divisible(a4 - a2*r, pi4)) continue;
            if (divisible(a6 + a4*r + a2*NTL::sqr(r) + NTL::power(r, 3), pi6)) {
                return true;
            }
        }
    }
    return false;
}

bool everywhere_minimal(const Coefficients& affine, const Coefficients& infinity,
                        unsigned q) {
    const bool cusp = q == 2
        ? NTL::IsZero(affine[0]) && NTL::IsZero(affine[2])
        : NTL::IsZero(affine[1]) && NTL::IsZero(affine[3]);
    if (cusp && (q == 2
        ? form_power(affine[3], 4) && form_power(affine[4], 2)
        : form_power(affine[4], 3))) return false;

    const zz_pX variable = polynomial({0, 1});
    if (nonminimal_at(infinity, variable, q)) return false;
    for (unsigned c = 0; c < q; ++c) {
        if (nonminimal_at(affine, polynomial({static_cast<int>(c), 1}), q)) return false;
    }
    // The degree bounds in the homogeneous local conditions exclude higher places.
    if (!cusp) {
        const std::vector<std::vector<int>> quadratics = q == 2
            ? std::vector<std::vector<int>>{{1, 1, 1}}
            : std::vector<std::vector<int>>{{1, 0, 1}, {2, 1, 1}, {2, 2, 1}};
        for (const auto& coefficients : quadratics) {
            if (nonminimal_at(affine, polynomial(coefficients), q)) return false;
        }
    }
    return true;
}

zz_pX quotient(const zz_pX& a, const zz_pX& pi, int n = 1) {
    zz_pX b = a, q, r;
    while (n--) {
        NTL::DivRem(q, r, b, pi);
        if (!NTL::IsZero(r)) throw std::runtime_error("inexact uniformizer division");
        b = q;
    }
    return b;
}

int valuation(zz_pX a, const zz_pX& pi) {
    if (NTL::IsZero(a)) return std::numeric_limits<int>::max();
    int v = 0;
    while (NTL::IsZero(a % pi)) {
        a = quotient(a, pi);
        ++v;
    }
    return v;
}

void transform(Coefficients& a, const zz_pX& r, const zz_pX& s, const zz_pX& t) {
    const auto [a1, a2, a3, a4, a6] = a;
    a[4] = a6 + r*(a4 + r*(a2 + r)) - t*(a3 + r*a1 + t);
    a[3] = a4 - s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r*r - 2*s*t;
    a[2] = a3 + r*a1 + 2*t;
    a[1] = a2 - s*a1 + 3*r - s*s;
    a[0] = a1 + 2*s;
}

struct ResidueField {
    using Element = NTL::zz_pE;
    int p, degree;
    long size;

    explicit ResidueField(const zz_pX& pi)
        : p(NTL::zz_p::modulus()), degree(NTL::deg(pi)), size(1) {
        NTL::zz_pE::init(pi);
        for (int i = 0; i < degree; ++i) size *= p;
    }

    Element reduce(const zz_pX& a) const { return NTL::conv<Element>(a); }
    zz_pX mod(const zz_pX& a) const { return NTL::rep(reduce(a)); }
    zz_pX inverse(const zz_pX& a) const { return NTL::rep(NTL::inv(reduce(a))); }
    zz_pX root(const zz_pX& a) const {
        return NTL::rep(NTL::power(reduce(a), size/p));
    }

    int character(const Element& a) const {
        const Element v = NTL::power(a, (size-1)/2);
        return NTL::IsOne(v) ? 1 : (NTL::IsZero(v) ? 0 : -1);
    }

    int trace(Element a) const {
        Element sum;
        for (int i = 0; i < degree; ++i) {
            sum += a;
            a = NTL::power(a, p);
        }
        return NTL::rep(NTL::coeff(NTL::rep(sum), 0));
    }

    bool quadratic(const Element& a, const Element& b, const Element& c) const {
        if (NTL::IsZero(a)) return !NTL::IsZero(b) || NTL::IsZero(c);
        if (p == 3) return character(b*b - 4*a*c) >= 0;
        if (NTL::IsZero(b)) return true;
        return trace(a*c/(b*b)) == 0;
    }

    bool quadratic(const zz_pX& a, const zz_pX& b, const zz_pX& c) const {
        return quadratic(reduce(a), reduce(b), reduce(c));
    }

    int cubic_roots(const zz_pX& b, const zz_pX& c, const zz_pX& d) const {
        NTL::zz_pEX f, x;
        NTL::SetCoeff(f, 3, 1);
        NTL::SetCoeff(f, 2, reduce(b));
        NTL::SetCoeff(f, 1, reduce(c));
        NTL::SetCoeff(f, 0, reduce(d));
        NTL::SetCoeff(x, 1, 1);
        return NTL::deg(NTL::GCD(f, NTL::PowerXMod(size, NTL::zz_pEXModulus(f)) - x));
    }
};

std::vector<int> component_cycles(const std::string& type, int n, int c, int trace,
                                  int degree) {
    std::vector<int> cycles;
    if (type == "III") cycles = {1};
    else if (type == "IV") cycles = c == 3 ? std::vector<int>{1, 1} : std::vector<int>{2};
    else if (type == "IV*") {
        cycles = c == 3 ? std::vector<int>(6, 1) : std::vector<int>{1, 1, 2, 2};
    } else if (type == "III*") cycles.assign(7, 1);
    else if (type == "II*") cycles.assign(8, 1);
    else if (type == "I" && n > 1) {
        if (trace == 1) cycles.assign(n-1, 1);
        else {
            cycles.assign((n-1)/2, 2);
            if (n % 2 == 0) cycles.push_back(1);
        }
    } else if (type == "I*") {
        if (n == 0 && c == 1) cycles = {1, 3};
        else if (c == 4) cycles.assign(n+4, 1);
        else {
            cycles.assign(n+2, 1);
            cycles.push_back(2);
        }
    }
    for (int& length : cycles) length *= degree;
    std::sort(cycles.begin(), cycles.end());
    return cycles;
}

LocalFibre local_reduction(Coefficients a, const zz_pX& pi, int vd,
                           const std::string& place) {
    NTL::zz_pEPush context;
    const ResidueField k(pi);
    using Element = ResidueField::Element;
    const int p = k.p;
    auto finish = [&](const std::string& type, int n, int components, int c, int trace = 0) {
        const int conductor = vd + 1 - components;
        if (conductor < 0) throw std::runtime_error("negative Ogg conductor");
        return LocalFibre{place, k.degree, type, n, vd, conductor, c, trace,
                          component_cycles(type, n, c, trace, k.degree)};
    };
    for (;;) {
        if (vd == 0) return finish("I", 0, 1, 1);
        WeierstrassInvariants b(a);
        const Element primary = k.reduce(p == 2 ? a[0] : b.b2);
        if (!NTL::IsZero(primary)) {
            const int trace = p == 3 ? k.character(primary)
                : (k.trace(k.reduce(a[1])/(primary*primary)
                    + k.reduce(a[2])/(primary*primary*primary)) ? -1 : 1);
            return finish("I", vd, vd, trace == 1 ? vd : (vd % 2 ? 1 : 2), trace);
        }
        zz_pX r, t;
        if (p == 2) {
            r = k.root(a[3]);
            t = k.root(((r + a[1])*r + a[3])*r + a[4]);
        } else {
            r = k.root(-b.b6);
            t = a[0]*r + a[2];
        }
        r = k.mod(r);
        t = k.mod(t);
        transform(a, r, zz_pX(), t);
        b = WeierstrassInvariants(a);
        if (valuation(a[4], pi) < 2) return finish("II", 0, 1, 1);
        if (valuation(b.b8, pi) < 3) return finish("III", 0, 2, 2);
        if (valuation(b.b6, pi) < 3) {
            return finish("IV", 0, 3,
                k.quadratic(zz_pX(1), quotient(a[2], pi), -quotient(a[4], pi, 2)) ? 3 : 1);
        }
        zz_pX s;
        if (p == 2) {
            s = k.root(a[1]);
            t = pi*k.root(quotient(a[4], pi, 2));
        } else {
            s = a[0];
            t = a[2];
        }
        transform(a, zz_pX(), s, t);
        const zz_pX cb = quotient(a[1], pi);
        const zz_pX cc = quotient(a[3], pi, 2);
        const zz_pX cd = quotient(a[4], pi, 3);
        const Element rb = k.reduce(cb), rc = k.reduce(cc), rd = k.reduce(cd);
        const Element w = 27*rd*rd - rb*rb*rc*rc + 4*rb*rb*rb*rd - 18*rb*rc*rd + 4*rc*rc*rc;
        const Element x = 3*rc - rb*rb;
        if (!NTL::IsZero(w)) return finish("I*", 0, 5, 1 + k.cubic_roots(cb, cc, cd));
        if (!NTL::IsZero(x)) {
            r = pi*(p == 2 ? k.root(cc) : k.mod(cc*k.inverse(cb)));
            transform(a, r, zz_pX(), zz_pX());
            int ix = 3, iy = 3;
            zz_pX mx = pi*pi, my = mx;
            for (;;) {
                zz_pX a2 = quotient(a[1], pi), a3 = quotient(a[2], my);
                zz_pX a4 = quotient(a[3], pi*mx), a6 = quotient(a[4], mx*my);
                const Element u = k.reduce(a3), v = k.reduce(a6);
                int c;
                if (NTL::IsZero(u*u + 4*v)) {
                    t = my*(p == 2 ? k.root(a6) : k.mod(a3));
                    transform(a, zz_pX(), zz_pX(), t);
                    my *= pi;
                    ++iy;
                    a2 = quotient(a[1], pi);
                    a4 = quotient(a[3], pi*mx);
                    a6 = quotient(a[4], mx*my);
                    const Element aa = k.reduce(a2), bb = k.reduce(a4), ccc = k.reduce(a6);
                    if (NTL::IsZero(bb*bb - 4*ccc*aa)) {
                        r = mx*(p == 2 ? k.root(a6*k.inverse(a2)) : k.mod(a4*k.inverse(a2)));
                        transform(a, r, zz_pX(), zz_pX());
                        mx *= pi;
                        ++ix;
                        if (ix + iy > vd + 5) {
                            throw std::runtime_error("Tate I* loop exceeded discriminant");
                        }
                        continue;
                    }
                    c = k.quadratic(aa, bb, ccc) ? 4 : 2;
                } else c = k.quadratic(Element(1), u, -v) ? 4 : 2;
                const int n = ix + iy - 5;
                return finish("I*", n, n+5, c);
            }
        }
        r = pi*(p == 2 ? k.mod(cb) : k.root(-cd));
        transform(a, r, zz_pX(), zz_pX());
        const zz_pX a3 = quotient(a[2], pi, 2), a6 = quotient(a[4], pi, 4);
        const Element u = k.reduce(a3), v = k.reduce(a6);
        if (!NTL::IsZero(u*u + 4*v)) {
            return finish("IV*", 0, 7, k.quadratic(Element(1), u, -v) ? 3 : 1);
        }
        t = pi*pi*(p == 2 ? -k.root(a6) : k.mod(a3));
        transform(a, zz_pX(), zz_pX(), t);
        if (valuation(a[3], pi) < 4) return finish("III*", 0, 8, 2);
        if (valuation(a[4], pi) < 6) return finish("II*", 0, 9, 1);
        constexpr int weights[] = {1, 2, 3, 4, 6};
        for (int j = 0; j < 5; ++j) a[j] = quotient(a[j], pi, weights[j]);
        vd -= 12;
        if (vd < 0) throw std::runtime_error("Tate scaling exceeded discriminant");
    }
}

Coefficients frobenius_twist(const Coefficients& a, unsigned p) {
    Coefficients result;
    for (int j = 0; j < 5; ++j) {
        for (long i = 0; i <= NTL::deg(a[j]); ++i) {
            NTL::SetCoeff(result[j], p*i, NTL::coeff(a[j], i));
        }
    }
    return result;
}

std::uint64_t place_code(const zz_pX& f, unsigned p) {
    std::uint64_t n = 0;
    for (long i = NTL::deg(f); i >= 0; --i) n = n*p + NTL::rep(NTL::coeff(f, i));
    return n;
}

void trim(std::vector<int>& coefficients) {
    while (coefficients.size() > 1 && coefficients.back() == 0) coefficients.pop_back();
}
}

Model model_from_key(unsigned q, std::uint64_t key) {
    if (q != 2 && q != 3) throw std::invalid_argument("q must be the prime 2 or 3");
    const std::uint64_t bound = q == 2 ? 4294967296ULL : 7625597484987ULL;
    if (key >= bound) throw std::invalid_argument("key is outside the encoding range");

    NTL::zz_pPush context(q);
    Model model;
    model.q = q;
    const std::array<bool, 5> included = q == 2
        ? std::array<bool, 5>{true, false, true, true, true}
        : std::array<bool, 5>{false, true, false, true, true};
    for (int i = 4; i >= 0; --i) {
        model.affine[i].resize(model_degrees[i] + 1);
        if (!included[i]) continue;
        for (int& digit : model.affine[i]) {
            digit = key % q;
            key /= q;
        }
    }
    model.key = 0;
    Coefficients affine;
    Coefficients infinity;
    for (int i = 0; i < 5; ++i) {
        model.infinity[i] = model.affine[i];
        std::reverse(model.infinity[i].begin(), model.infinity[i].end());
        if (included[i]) {
            for (int digit : model.infinity[i]) model.key = q*model.key + digit;
        }
        affine[i] = polynomial(model.affine[i]);
        infinity[i] = polynomial(model.infinity[i]);
        trim(model.affine[i]);
        trim(model.infinity[i]);
    }
    model.kind = fibre_kind(affine, q);
    model.minimality = everywhere_minimal(affine, infinity, q)
        ? "everywhere minimal" : "nonminimal somewhere";
    model.status = model.kind == "smooth" && model.minimality == "everywhere minimal"
        ? "partial" : "unsupported";
    return model;
}

LocalModels local_models(const Model& model) {
    if (model.q != 2 && model.q != 3) throw std::invalid_argument("q must be the prime 2 or 3");
    if (model.kind != "smooth" || model.minimality != "everywhere minimal") {
        throw std::invalid_argument("local arithmetic requires an everywhere minimal elliptic K3");
    }
    NTL::zz_pPush context(model.q);
    Coefficients affine, infinity;
    for (int j = 0; j < 5; ++j) {
        affine[j] = polynomial(model.affine[j]);
        infinity[j] = polynomial(model.infinity[j]);
    }
    zz_pX delta = WeierstrassInvariants(affine).discriminant();
    if (NTL::IsZero(delta)) throw std::runtime_error("singular generic fibre");
    const int infinity_valuation = 24 - NTL::deg(delta);
    if (infinity_valuation < 0) throw std::runtime_error("discriminant exceeds K3 chart degree");
    const auto twisted = frobenius_twist(affine, model.q);
    NTL::MakeMonic(delta);
    NTL::vec_pair_zz_pX_long factors;
    NTL::CanZass(factors, delta);
    LocalModels result;
    auto append = [&](const Coefficients& a, const Coefficients& twist, const zz_pX& pi,
                      int vd, const std::string& place) {
        try {
            const auto original = local_reduction(a, pi, vd, place);
            if (original.discriminant_valuation != vd) {
                throw std::runtime_error("input model is not locally minimal");
            }
            const auto other = local_reduction(twist, pi, model.q*vd, place);
            if (original.conductor) result.fibres.push_back(original);
            if (other.conductor) result.twist_fibres.push_back(other);
        } catch (const std::exception& error) {
            throw std::runtime_error("local place " + place + ": " + error.what());
        }
    };
    for (long i = 0; i < factors.length(); ++i) {
        const auto& pi = factors[i].a;
        append(affine, twisted, pi, factors[i].b, std::to_string(place_code(pi, model.q)));
    }
    if (infinity_valuation > 0) {
        // Twist the weight-2i chart before minimization; its weights become 2ip.
        append(infinity, frobenius_twist(infinity, model.q), polynomial({0, 1}),
               infinity_valuation, "inf");
    }
    auto order = [](const LocalFibre& a, const LocalFibre& b) {
        if (a.place == "inf" || b.place == "inf") return b.place == "inf" && a.place != "inf";
        return std::stoull(a.place) < std::stoull(b.place);
    };
    std::sort(result.fibres.begin(), result.fibres.end(), order);
    std::sort(result.twist_fibres.begin(), result.twist_fibres.end(), order);
    return result;
}

}
