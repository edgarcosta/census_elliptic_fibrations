#include "invariants.hpp"
#include "../records/format.hpp"
#include "../models/weierstrass.hpp"

#include <NTL/lzz_pX.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_factor.h>
#include <gmpxx.h>

#include <algorithm>
#include <numeric>
#include <ostream>
#include <set>
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
    NTL::zz_pPush context(model.q);
    Q maximum = 0;
    for (const auto& fibre : local.fibres) {
        maximum += fibre.degree*maximum_correction(fibre);
    }
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
        result.points.push_back(p.identity ? "O" : polynomial_code(p.x.n, model.q)+":"+polynomial_code(p.x.d, model.q)
            +":"+polynomial_code(p.y.n, model.q)+":"+polynomial_code(p.y.d, model.q));
    }
    const int first = points.size()/exponent;
    if (points.size() % exponent || exponent % first) throw std::runtime_error("invalid torsion invariant factors");
    if (first > 1) result.factors.push_back(first);
    if (exponent > 1) result.factors.push_back(exponent);
    return result;
}

namespace {
using Z = mpz_class;
using Poly = std::vector<Z>;

struct IntegerPolynomial {
    fmpz_poly_t value;

    explicit IntegerPolynomial(const Poly& coefficients = {}) {
        fmpz_poly_init2(value, coefficients.size());
        for (std::size_t i = 0; i < coefficients.size(); ++i) fmpz_set_mpz(value->coeffs+i, coefficients[i].get_mpz_t());
        _fmpz_poly_set_length(value, coefficients.size());
        _fmpz_poly_normalise(value);
    }
    ~IntegerPolynomial() { fmpz_poly_clear(value); }
    IntegerPolynomial(const IntegerPolynomial&) = delete;
    IntegerPolynomial& operator=(const IntegerPolynomial&) = delete;

    Poly coefficients() const {
        Poly result(fmpz_poly_length(value));
        for (std::size_t i = 0; i < result.size(); ++i) {
            fmpz_get_mpz(result[i].get_mpz_t(), value->coeffs+i);
        }
        return result;
    }
};

Z integer_power(const Z& base, int exponent) {
    Z value;
    mpz_pow_ui(value.get_mpz_t(), base.get_mpz_t(), exponent);
    return value;
}

void multiply_binomial(Poly& a, int n, const Z& q) {
    const Z scale = integer_power(q, n);
    const std::size_t size = a.size();
    a.resize(size+n);
    for (std::size_t i = size; i-- > 0;) a[i+n] -= scale*a[i];
}

Q evaluate_inverse(const IntegerPolynomial& polynomial, const Z& q) {
    Q result = 0;
    const auto coefficients = polynomial.coefficients();
    for (auto i = coefficients.rbegin(); i != coefficients.rend(); ++i) result = result/q+*i;
    return result;
}

void include_order(Z& degree, int n) {
    mpz_lcm_ui(degree.get_mpz_t(), degree.get_mpz_t(), n);
}

const std::vector<std::pair<int, Poly>>& cyclotomic_factors() {
    static const auto factors = [] {
        std::vector<std::pair<int, Poly>> result;
        // phi(n) >= sqrt(n/2); a factor of an L-polynomial has degree <= 20.
        fmpz_poly_t f;
        fmpz_poly_init(f);
        for (int n = 1; n <= 2*20*20; ++n) {
            fmpz_poly_cyclotomic(f, n);
            const int degree = fmpz_poly_degree(f);
            if (degree > 20) continue;
            Poly c(degree+1);
            for (int i = 0; i <= degree; ++i) c[i] = fmpz_poly_get_coeff_si(f, i);
            if (c[0] < 0) for (auto& value : c) value = -value;
            result.emplace_back(n, std::move(c));
        }
        fmpz_poly_clear(f);
        return result;
    }();
    return factors;
}
}

DerivedInvariants derived_invariants(const Model& model, const LocalModels& local,
    const TorsionData& torsion, const LPolynomial& lfun) {
    const Poly& L = lfun.coefficients;
    if (L.empty() || L[0] != 1 || L.size() > 21) throw std::invalid_argument("invalid K3 L-polynomial");
    const Z q = model.q;
    DerivedInvariants result;
    result.p2 = L;
    multiply_binomial(result.p2, 1, q);
    multiply_binomial(result.p2, 1, q);
    int cycle_count = 0, cycle_sum = 0;
    Z tamagawa = 1, torsion_order = 1;
    for (int n : torsion.factors) torsion_order *= n;
    for (const auto& f : local.fibres) {
        tamagawa *= f.tamagawa;
        for (int n : f.component_cycles) {
            multiply_binomial(result.p2, n, q);
            ++cycle_count;
            cycle_sum += n;
            include_order(result.ns_splitting_degree, n);
        }
    }
    if (result.p2.size() != 23) throw std::invalid_argument("L and component cycles must give P2 of degree 22");
    const IntegerPolynomial linear(Poly{1, -q});
    IntegerPolynomial quotient(L), divided;
    while (fmpz_poly_divides(divided.value, quotient.value, linear.value)) {
        fmpz_poly_swap(quotient.value, divided.value);
        ++result.rank;
    }
    result.picard_number = 2+result.rank+cycle_count;
    IntegerPolynomial remainder(L);
    for (const auto& [order, coefficients] : cyclotomic_factors()) {
        if (slong(coefficients.size()) > fmpz_poly_length(remainder.value)) continue;
        Poly scaled_coefficients = coefficients;
        Z scale = 1;
        for (auto& c : scaled_coefficients) { c *= scale; scale *= q; }
        const IntegerPolynomial scaled(scaled_coefficients);
        while (fmpz_poly_divides(divided.value, remainder.value, scaled.value)) {
            fmpz_poly_swap(remainder.value, divided.value);
            result.geometric_rank += int(coefficients.size())-1;
            include_order(result.ns_splitting_degree, order);
        }
    }
    result.geometric_picard_number = 2+result.geometric_rank+cycle_sum;
    Q slope = 1;
    for (std::size_t i = 1; i < L.size(); ++i) {
        if (L[i] == 0) continue;
        Z c = L[i];
        int v = 0;
        while (c % model.q == 0) { c /= model.q; ++v; }
        slope = std::min(slope, Q(Q(v)/i));
    }
    if (slope < 1) {
        const Q height = 1/(1-slope);
        if (height.get_den() != 1 || height < 1 || height > 10) throw std::runtime_error("invalid formal Brauer height");
        result.brauer_height = height.get_num().get_si();
    }
    IntegerPolynomial algebraic_quotient(result.p2);
    for (int i = 0; i < result.picard_number; ++i) {
        if (!fmpz_poly_divides(divided.value, algebraic_quotient.value, linear.value)) throw std::runtime_error("Picard multiplicity mismatch");
        fmpz_poly_swap(algebraic_quotient.value, divided.value);
    }
    const Q artin_tate = q*evaluate_inverse(algebraic_quotient, q);
    if (artin_tate.get_den() != 1 || artin_tate <= 0) throw std::runtime_error("invalid Artin-Tate product");
    result.artin_tate_product = artin_tate.get_num();
    result.sha_reg = q*evaluate_inverse(quotient, q)*torsion_order*torsion_order/tamagawa;
    const Z leading_scale = integer_power(q, 22);
    if (result.p2.back() != leading_scale && result.p2.back() != -leading_scale) {
        throw std::runtime_error("invalid cohomology determinant");
    }
    result.zeta_sign = result.p2.back() == leading_scale ? 1 : -1;
    IntegerPolynomial characteristic(result.p2), power_sums;
    fmpz_poly_reverse(characteristic.value, characteristic.value, result.p2.size());
    fmpz_poly_power_sums(power_sums.value, characteristic.value, 6);
    const Poly traces = power_sums.coefficients();
    for (int n = 1; n <= 5; ++n) {
        result.point_counts.push_back(1+integer_power(q, 2*n)+(n < int(traces.size()) ? traces[n] : Z(0)));
    }
    NTL::zz_pPush context(model.q);
    const WeierstrassInvariants w(polynomials(model.affine));
    P numerator = NTL::power(w.b2*w.b2-24*w.b4, 3), denominator = w.discriminant();
    const P common = NTL::GCD(numerator, denominator);
    numerator /= common;
    denominator /= common;
    result.j_degree = std::max(NTL::deg(numerator), NTL::deg(denominator));
    if (result.j_degree > 0) {
        int common_exponent = 0;
        for (const auto& f : {numerator, denominator}) {
            for (int i = 1; i <= NTL::deg(f); ++i) {
                if (!NTL::IsZero(NTL::coeff(f, i))) common_exponent = std::gcd(common_exponent, i);
            }
        }
        while (common_exponent % model.q == 0) {
            ++result.j_inseparability_exponent;
            common_exponent /= model.q;
        }
    }
    return result;
}

void write_derived_invariants(std::ostream& output, const DerivedInvariants& r) {
    output << r.rank << ',' << r.geometric_rank << ',' << r.picard_number
           << ',' << r.geometric_picard_number << ',';
    write_joined(output, r.p2, ' ');
    output << ',' << r.zeta_sign << ',';
    write_joined(output, r.point_counts, ' ');
    output << ',';
    if (r.brauer_height) output << r.brauer_height;
    else output << "infinity";
    output << ',' << r.artin_tate_product << ',' << r.ns_splitting_degree
           << ',' << r.j_degree << ',' << r.j_inseparability_exponent << ','
           << r.sha_reg.get_num();
    if (r.sha_reg.get_den() != 1) output << '/' << r.sha_reg.get_den();
}

namespace {
std::vector<Z> square_divisors(const Z& value) {
    fmpz_t n;
    fmpz_init(n);
    fmpz_set_mpz(n, value.get_mpz_t());
    fmpz_factor_t factors;
    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    std::vector<Z> squares{1};
    for (slong i = 0; i < factors->num; ++i) {
        Z prime;
        fmpz_get_mpz(prime.get_mpz_t(), factors->p+i);
        const std::size_t count = squares.size();
        Z power = 1;
        for (ulong j = 1; j <= factors->exp[i]/2; ++j) {
            power *= prime*prime;
            for (std::size_t k = 0; k < count; ++k) squares.push_back(squares[k]*power);
        }
    }
    fmpz_factor_clear(factors);
    fmpz_clear(n);
    std::sort(squares.begin(), squares.end());
    return squares;
}

std::set<Q> rational_corrections(const LocalFibre& f) {
    std::set<Q> values{0};
    if (f.type == "I") {
        if (f.tamagawa == f.index) {
            for (int i = 0; i < f.index; ++i) values.insert(fraction(i*(f.index-i), f.index));
        } else if (f.index % 2 == 0) values.insert(fraction(f.index, 4));
    } else if (f.type == "I*") {
        if (f.tamagawa == 2 || f.tamagawa == 4) values.insert(1);
        if (f.tamagawa == 4) values.insert(1+fraction(f.index, 4));
    } else if (f.type == "III") values.insert(fraction(1, 2));
    else if (f.type == "III*") values.insert(fraction(3, 2));
    else if (f.type == "IV" && f.tamagawa == 3) values.insert(fraction(2, 3));
    else if (f.type == "IV*" && f.tamagawa == 3) values.insert(fraction(4, 3));
    std::set<Q> weighted;
    for (const auto& c : values) weighted.insert(f.degree*c);
    return weighted;
}

Z rational_floor(const Q& value) {
    Z result;
    mpz_fdiv_q(result.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());
    return result;
}

Q rational_power(const Q& value, int exponent) {
    return Q(integer_power(value.get_num(), exponent))/integer_power(value.get_den(), exponent);
}
}

CandidateInvariants admissible_candidates(const DerivedInvariants& data,
    const std::vector<LocalFibre>& fibres) {
    if (data.rank < 0 || data.rank > 20 || data.picard_number < 2 || data.picard_number > 22
        || data.artin_tate_product <= 0 || data.sha_reg <= 0) {
        throw std::invalid_argument("candidates require K3 ranks and positive exact products");
    }
    std::set<Q> corrections{0};
    if (data.rank > 0) {
        for (const auto& f : fibres) {
            std::set<Q> sums;
            const auto local = rational_corrections(f);
            for (const auto& a : corrections) for (const auto& b : local) sums.insert(a+b);
            corrections = std::move(sums);
        }
    }
    Q lower_bound = 0;
    if (data.rank > 1) {
        Q minimum = 4;
        for (const auto& c : corrections) {
            const Z intersections = std::max(Z(0), Z(rational_floor((c-4)/2)+1));
            minimum = std::min(minimum, Q(4+2*intersections-c));
        }
        const Q gamma_powers[] = {1, fraction(4, 3), 2, 4, 8};
        const Q gamma = data.rank <= 5 ? gamma_powers[data.rank-1]
            : rational_power(fraction(4, 3), data.rank*(data.rank-1)/2);
        lower_bound = rational_power(minimum, data.rank)/gamma;
    }
    CandidateInvariants result;
    for (const auto& s : square_divisors(data.artin_tate_product)) {
        const Z discriminant = data.artin_tate_product/s;
        const bool parity = data.picard_number % 2 ? discriminant % 2 == 0
            : discriminant % 4 == 0 || discriminant % 4 == (data.picard_number % 4 == 0 ? 3 : 1);
        if (!parity) continue;
        const Q regulator = data.sha_reg/s;
        if (data.rank == 0) {
            if (regulator != 1) continue;
        } else if (data.rank == 1) {
            const bool height = std::any_of(corrections.begin(), corrections.end(), [&](const Q& c) {
                const Q intersections = (regulator+c-4)/2;
                return intersections >= 0 && intersections.get_den() == 1;
            });
            if (!height) continue;
        } else if (regulator < lower_bound) continue;
        result.sha.push_back(s);
        result.regulators.push_back(regulator);
    }
    return result;
}

void write_candidate_invariants(std::ostream& output, const CandidateInvariants& candidates) {
    write_list(output, candidates.sha);
    output << ',';
    write_list(output, candidates.regulators);
}

}
