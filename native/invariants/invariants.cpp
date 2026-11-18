#include "invariants.hpp"
#include "../models/weierstrass.hpp"

#include <NTL/lzz_pX.h>
#include <gmpxx.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace census {
namespace {
using P = NTL::zz_pX;
using Q = mpq_class;

Q fraction(int n, int d) { return Q(n)/d; }

P decode(std::uint64_t n, unsigned q) {
    P f;
    for (int i = 0; n; ++i, n /= q) NTL::SetCoeff(f, i, n % q);
    return f;
}

struct Rational {
    P n, d;
    Rational(long value = 0) : n(value), d(1) {}
    Rational(const P& value) : n(value), d(1) {}
    Rational(P numerator, P denominator) : n(std::move(numerator)), d(std::move(denominator)) {
        if (NTL::IsZero(d)) throw std::runtime_error("zero rational denominator");
        const P common = NTL::GCD(n, d);
        n /= common;
        d /= common;
        const auto leading = NTL::LeadCoeff(d);
        n /= leading;
        d /= leading;
    }
};

Rational operator+(const Rational& a, const Rational& b) { return {a.n*b.d+b.n*a.d, a.d*b.d}; }
Rational operator-(const Rational& a, const Rational& b) { return {a.n*b.d-b.n*a.d, a.d*b.d}; }
Rational operator-(const Rational& a) { return {-a.n, a.d}; }
Rational operator*(const Rational& a, const Rational& b) { return {a.n*b.n, a.d*b.d}; }
Rational operator/(const Rational& a, const Rational& b) { return {a.n*b.d, a.d*b.n}; }
bool operator==(const Rational& a, const Rational& b) { return a.n == b.n && a.d == b.d; }

struct Point {
    Rational x, y;
    bool identity = true;
    Point() = default;
    Point(Rational x0, Rational y0) : x(std::move(x0)), y(std::move(y0)), identity(false) {}
};

bool operator==(const Point& a, const Point& b) {
    if (a.identity || b.identity) return a.identity == b.identity;
    return a.x == b.x && a.y == b.y;
}

using census::valuation;

int valuation(const Rational& f, const P& pi) {
    return NTL::IsZero(f.n) ? infinite_valuation : valuation(f.n, pi)-valuation(f.d, pi);
}

Rational reciprocal(const Rational& f, int weight) {
    if (NTL::IsZero(f.n)) return f;
    P n = NTL::reverse(f.n), d = NTL::reverse(f.d);
    const int shift = weight + NTL::deg(f.d) - NTL::deg(f.n);
    if (shift >= 0) NTL::LeftShift(n, n, shift);
    else NTL::LeftShift(d, d, -shift);
    return {n, d};
}

struct HeightChart : WeierstrassInvariants {
    std::array<Rational, 5> a;
    P c4;

    explicit HeightChart(const Coefficients& coefficients)
        : WeierstrassInvariants(coefficients), c4(b2*b2-24*b4) {
        for (int i = 0; i < 5; ++i) a[i] = Rational(coefficients[i]);
    }

    Q correction(const Rational& x, const Rational& y, const P& pi, int vd) const {
        const auto& [a1, a2, a3, a4, a6] = a;
        const int f2 = valuation(Rational(2)*y+a1*x+a3, pi);
        if (valuation(Rational(3)*x*x+Rational(2)*a2*x+a4-a1*y, pi) <= 0 || f2 <= 0) return 0;
        if (valuation(c4, pi) == 0) {
            const Q alpha = std::min(fraction(1, 2), fraction(f2, vd));
            return vd*alpha*(1-alpha);
        }
        const int f3 = valuation(Rational(3)*x*x*x*x+b2*x*x*x
            +Rational(3)*b4*x*x+Rational(3)*b6*x+b8, pi);
        return f3 >= 3*f2 ? fraction(2*f2, 3) : fraction(f3, 4);
    }
};

struct Curve {
    HeightChart affine, infinity;
    struct Place { P pi; int degree, vd; bool infinity; };
    std::vector<Place> places;

    Curve(const Model& model, const LocalModels& local)
        : affine(polynomials(model.affine)), infinity(polynomials(model.infinity)) {
        for (const auto& row : local.fibres) {
            const bool inf = row.place == "inf";
            places.push_back({inf ? polynomial({0, 1}) : decode(std::stoull(row.place), model.q),
                              row.degree, row.discriminant_valuation, inf});
        }
    }

    Point negate(const Point& p) const {
        return p.identity ? p : Point(p.x, -p.y-affine.a[0]*p.x-affine.a[2]);
    }

    Point add(const Point& p, const Point& s) const {
        if (p.identity) return s;
        if (s.identity) return p;
        const auto& [a1, a2, a3, a4, a6] = affine.a;
        Rational slope;
        if (p.x == s.x) {
            if (s == negate(p)) return {};
            slope = (Rational(3)*p.x*p.x+Rational(2)*a2*p.x+a4-a1*p.y)
                /(Rational(2)*p.y+a1*p.x+a3);
        } else slope = (s.y-p.y)/(s.x-p.x);
        const Rational intercept = p.y-slope*p.x;
        const Rational x = slope*slope+a1*slope-a2-p.x-s.x;
        return {x, -(slope+a1)*x-intercept-a3};
    }

    Q height(const Point& p) const {
        if (p.identity) return 0;
        Q value = 4+NTL::deg(p.x.d)+std::max(0L, NTL::deg(p.x.n)-NTL::deg(p.x.d)-4);
        for (const auto& place : places) {
            const Rational x = place.infinity ? reciprocal(p.x, 4) : p.x;
            const Rational y = place.infinity ? reciprocal(p.y, 6) : p.y;
            value -= place.degree*(place.infinity ? infinity : affine).correction(x, y, place.pi, place.vd);
        }
        if (value < 0) throw std::runtime_error("negative section height");
        return value;
    }
};

Q maximum_correction(const LocalFibre& f) {
    if (f.type == "I") return fraction((f.index/2)*(f.index-f.index/2), f.index);
    if (f.type == "I*") return fraction(f.index+4, 4);
    if (f.type == "III") return fraction(1, 2);
    if (f.type == "IV") return fraction(2, 3);
    if (f.type == "IV*") return fraction(4, 3);
    if (f.type == "III*") return fraction(3, 2);
    if (f.type == "II" || f.type == "II*") return 0;
    throw std::runtime_error("unknown Kodaira type in torsion bound");
}

std::vector<P> quadratic_roots(const P& linear, const P& rhs, unsigned q, int degree) {
    if (q == 3) {
        const P square = linear*linear+rhs;
        if (NTL::IsZero(square)) return {linear};
        const int n = NTL::deg(square)/2;
        if (NTL::deg(square) % 2 || NTL::LeadCoeff(square) != 1) return {};
        P root;
        NTL::SetCoeff(root, n, 1);
        for (int j = n-1; j >= 0; --j) {
            NTL::SetCoeff(root, j, -NTL::coeff(square-root*root, n+j));
        }
        if (root*root != square) return {};
        return {linear+root, linear-root};
    }
    // Over GF(2), B -> B^2+linear*B is linear, with kernel {0, linear}.
    std::vector<P> basis(2*degree+1), preimages(2*degree+1);
    for (int j = 0; j <= degree; ++j) {
        P preimage;
        NTL::SetCoeff(preimage, j, 1);
        P image = preimage*preimage+linear*preimage;
        while (!NTL::IsZero(image)) {
            const int k = NTL::deg(image);
            if (NTL::IsZero(basis[k])) {
                basis[k] = image;
                preimages[k] = preimage;
                break;
            }
            image += basis[k];
            preimage += preimages[k];
        }
    }
    P remainder = rhs, root;
    while (!NTL::IsZero(remainder)) {
        const int k = NTL::deg(remainder);
        if (k >= int(basis.size()) || NTL::IsZero(basis[k])) return {};
        remainder += basis[k];
        root += preimages[k];
    }
    if (NTL::IsZero(linear)) return {root};
    return {root, root+linear};
}

std::vector<Point> zero_height_sections(const Curve& curve, unsigned q, int bound) {
    std::vector<Point> points(1);
    if (bound < 0) return points;
    std::vector<P> denominators{P(1)};
    if (bound == 1) {
        for (unsigned c = 0; c < q; ++c) denominators.push_back(polynomial({int(c), 1}));
    }
    const auto& a = curve.affine.a;
    unsigned count = 1;
    for (int i = 0; i <= 4+2*bound; ++i) count *= q;
    for (const P& g : denominators) {
        const P g2 = g*g, g3 = g2*g;
        for (unsigned code = 0; code < count; ++code) {
            const P A = decode(code, q);
            if (!NTL::IsOne(NTL::GCD(A, g))) continue;
            // Poles at infinity count towards (P.O), including polynomial x.
            const int excess = NTL::deg(A)-2*NTL::deg(g)-4;
            if (NTL::deg(g)+std::max(0, (excess+1)/2) > bound) continue;
            const P rhs = A*A*A+a[1].n*A*A*g2+a[3].n*A*g2*g2+a[4].n*g3*g3;
            const P linear = a[0].n*A*g+a[2].n*g3;
            for (const P& B : quadratic_roots(linear, rhs, q, 6+3*bound)) {
                if (NTL::deg(B) > 6+3*bound) continue;
                const Point point(Rational(A, g2), Rational(B, g3));
                if (curve.height(point) == 0) points.push_back(point);
            }
        }
    }
    return points;
}
}

TorsionData complete_torsion(const Model& model, const LocalModels& local) {
    if (model.status != "partial") throw std::invalid_argument("torsion requires an everywhere minimal elliptic K3");
    NTL::zz_pPush context(model.q);
    Q maximum = 0;
    int discriminant_degree = 0;
    for (const auto& fibre : local.fibres) {
        maximum += fibre.degree*maximum_correction(fibre);
        discriminant_degree += fibre.degree*fibre.discriminant_valuation;
    }
    if (discriminant_degree != 24 || maximum > 6) throw std::runtime_error("K3 torsion bound violated");
    const int bound = maximum < 4 ? -1 : maximum < 6 ? 0 : 1;
    const Curve curve(model, local);
    const auto points = zero_height_sections(curve, model.q, bound);
    TorsionData result;
    int exponent = 1;
    for (const auto& p : points) {
        Point sum;
        int order = 0;
        do {
            sum = curve.add(sum, p);
            ++order;
            if (order > int(points.size()) || std::find(points.begin(), points.end(), sum) == points.end()) {
                throw std::runtime_error("torsion enumeration is not closed");
            }
        } while (!sum.identity);
        exponent = std::lcm(exponent, order);
        if (exponent > int(points.size())) throw std::runtime_error("torsion exponent exceeds group order");
        result.orders.push_back(order);
        result.points.push_back(p.identity ? "O" : polynomial_code(p.x.n, model.q)+":"+polynomial_code(p.x.d, model.q)
            +":"+polynomial_code(p.y.n, model.q)+":"+polynomial_code(p.y.d, model.q));
    }
    const int first = points.size()/exponent;
    if (points.size() % exponent || exponent % first) throw std::runtime_error("invalid torsion invariant factors");
    if (first > 1) result.factors.push_back(first);
    if (exponent > 1) result.factors.push_back(exponent);
    return result;
}

}
