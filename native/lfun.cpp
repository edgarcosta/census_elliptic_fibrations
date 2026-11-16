#include "lfun.hpp"
#include "weierstrass.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#if !defined(CENSUS_SCALAR_ONLY) && (defined(__x86_64__) || defined(__i386__)) \
    && (defined(__GNUC__) || defined(__clang__))
#define CENSUS_AVX2 1
#include <immintrin.h>
#endif

namespace census {

SignData functional_equation_sign(const LocalModels& local, unsigned q) {
    if (q != 2 && q != 3) throw std::invalid_argument("q must be the prime 2 or 3");
    auto summarize = [q](const std::vector<LocalFibre>& fibres) {
        std::array<int, 3> values{};
        for (const auto& fibre : fibres) {
            values[0] += fibre.degree*fibre.discriminant_valuation;
            values[1] += fibre.degree*fibre.conductor;
            if (fibre.tamagawa <= 0) throw std::runtime_error("nonpositive Tamagawa number");
            // Each closed place contributes once to C, without a degree exponent.
            for (int c = fibre.tamagawa; c % q == 0; c /= q) ++values[2];
        }
        if (values[0] % 12) throw std::runtime_error("nonintegral elliptic surface height");
        return values;
    };
    const auto original = summarize(local.fibres);
    const auto twist = summarize(local.twist_fibres);
    if (original[0] != 24) throw std::runtime_error("K3 discriminant degree must be 24");
    SignData result;
    result.degree = original[1] - 4;
    if (result.degree < 0 || result.degree > 20) {
        throw std::runtime_error("L-polynomial degree outside K3 range");
    }
    result.height = original[0]/12;
    result.twist_height = twist[0]/12;
    result.tamagawa_exponent = twist[2] - original[2];
    const int exponent = result.twist_height - result.height + result.tamagawa_exponent;
    result.sign = exponent % 2 ? -1 : 1;
    return result;
}

namespace {
constexpr std::int64_t modulus = 1062907;
constexpr std::int64_t omega = 662797;

std::int64_t modular_power(std::int64_t a, int e) {
    std::int64_t result = 1;
    for (; e; e >>= 1, a = a * a % modulus) {
        if (e & 1) result = result * a % modulus;
    }
    return result;
}

void radix3(std::vector<std::int64_t>& values, std::int64_t root) {
    const auto root2 = root * root % modulus;
    const int size = values.size();
    for (int stride = 1; stride < size; stride *= 3) {
        for (int block = 0; block < size; block += 3 * stride) {
            for (int j = block; j < block + stride; ++j) {
                const auto x = values[j], y = values[j + stride], z = values[j + 2 * stride];
                values[j] = (x + y + z) % modulus;
                values[j + stride] = (x + root * y + root2 * z) % modulus;
                values[j + 2 * stride] = (x + root2 * y + root * z) % modulus;
            }
        }
    }
}

void ordinary3(TraceTable& table) {
    const auto& field = table.field;
    bool prime = true;
    for (int divisor = 2; std::int64_t(divisor) * divisor <= modulus; ++divisor) {
        if (modulus % divisor == 0) prime = false;
    }
    if (!prime || modulus % 3 != 1 || modulus <= 2 * field.size || omega == 1
        || modular_power(omega, 3) != 1 || modular_power(omega, 2) != 400109) {
        throw std::runtime_error("invalid trace-table certificate constants");
    }
    std::vector<std::int64_t> h(field.size), chi(field.size);
    for (int e = 0; e < field.order; ++e) {
        const auto b = field.exponents[e];
        ++h[field.index(field.negate(field.add(field.exp(3LL * e), field.exp(2LL * e))))];
        chi[field.index(b)] = e & 1 ? modulus - 1 : 1;
    }
    ++h[0];
    std::int64_t sum = 0;
    for (const auto value : h) sum += value;
    if (h[0] != 2 || sum != field.size) throw std::runtime_error("invalid convolution weights");
    radix3(h, omega);
    radix3(chi, omega);
    for (int i = 0; i < field.size; ++i) h[i] = h[i] * chi[i] % modulus;
    radix3(h, omega * omega % modulus);
    const auto inverse = modular_power(field.size, modulus - 2);
    if (inverse * field.size % modulus != 1) throw std::runtime_error("invalid normalization");
    // M > 2Q makes centered recovery unique in [-Q,Q]; see native-mathematics.md.
    for (auto& value : h) {
        value = value * inverse % modulus;
        if (value > (modulus - 1) / 2) value -= modulus;
        if (value < -field.size || value > field.size) {
            throw std::runtime_error("trace-table certificate bound failed");
        }
    }
    table.ordinary.resize(field.order);
    for (int e = 0; e < field.order; ++e) table.ordinary[e] = -h[field.index(field.exponents[e])];
    table.nodal = -h[0];
}

void ordinary2(TraceTable& table) {
    const auto& field = table.field;
    std::vector<int> h(field.size);
    for (int e = 0; e < field.order; ++e) {
        h[field.exp(-2LL * e)] = field.trace(field.exponents[e]) ? -1 : 1;
    }
    for (int stride = 1; stride < field.size; stride *= 2) {
        for (int block = 0; block < field.size; block += 2 * stride) {
            for (int j = block; j < block + stride; ++j) {
                const int x = h[j], y = h[j + stride];
                h[j] = x + y;
                h[j + stride] = x - y;
            }
        }
    }
    table.ordinary.resize(field.order);
    for (int e = 0; e < field.order; ++e) {
        int dual = 0;
        for (int j = 0; j < field.degree; ++j) {
            dual |= field.trace(field.multiply(field.exponents[e], 1u << j)) << j;
        }
        table.ordinary[e] = -h[dual];
    }
    table.nodal = -h[0];
}

std::string canonical_table(const TraceTable& table) {
    std::ostringstream output;
    output << table.field.p << ' ' << table.field.degree << '\n';
    const auto line = [&output](const auto& values) {
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i) output << ' ';
            output << values[i];
        }
        output << '\n';
    };
    line(table.field.polynomial);
    line(table.ordinary);
    return output.str();
}

void check_cache(const TraceTable& table, const std::filesystem::path& directory) {
    const auto path = directory / ("f" + std::to_string(table.field.p) + "-d"
        + std::to_string(table.field.degree) + ".txt");
    const auto expected = canonical_table(table);
    if (std::filesystem::exists(path)) {
        std::ifstream input(path, std::ios::binary);
        if (!input) throw std::runtime_error("could not read trace-table cache: " + path.string());
        const std::string cached{std::istreambuf_iterator<char>(input), {}};
        if (input.bad()) throw std::runtime_error("could not read trace-table cache: " + path.string());
        if (cached != expected) {
            throw std::runtime_error("trace-table certificate mismatch: " + path.string());
        }
    } else {
        std::filesystem::create_directories(directory);
        std::ofstream output(path, std::ios::binary);
        output << expected;
        output.close();
        if (!output) throw std::runtime_error("could not write trace-table cache: " + path.string());
    }
}
}

TraceTable::TraceTable(int p, int degree, const std::filesystem::path& cache)
    : field(p, degree), nodal(0) {
    if (p == 3) ordinary3(*this);
    else ordinary2(*this);
    if (nodal != 1) throw std::runtime_error("nodal trace is not one");
    if (!cache.empty()) check_cache(*this, cache);
}

namespace {
struct CountingTable {
    TraceTable traces;
    std::array<std::vector<std::uint32_t>, 25> powers;
    std::array<int, 3> gauss{};

    CountingTable(int p, int degree) : traces(p, degree) {
        const auto& f = traces.field;
        std::vector<std::uint32_t> representatives;
        for (int e = 0; e < f.order; ++e) {
            int next = e, length = 0, smallest = e;
            do {
                smallest = std::min(smallest, next);
                next = p*next % f.order;
                ++length;
            } while (next != e);
            if (length == degree && smallest == e) representatives.push_back(f.exponents[e]);
            if (p == 3) gauss[f.trace(f.exponents[e])] += e & 1 ? -1 : 1;
        }
        if (degree == 1) representatives.push_back(0);
        for (auto& column : powers) column.resize(representatives.size());
        for (std::size_t i = 0; i < representatives.size(); ++i) {
            std::uint32_t value = 1;
            for (auto& column : powers) {
                column[i] = value;
                value = f.multiply(value, representatives[i]);
            }
        }
    }
};

void evaluate_scalar(const CountingTable& table, const std::vector<int>& poly,
                     std::vector<std::uint32_t>& output) {
    output.assign(table.powers[0].size(), 0);
    const auto& f = table.traces.field;
    for (std::size_t j = 0; j < poly.size(); ++j) {
        if (!poly[j]) continue;
        for (std::size_t i = 0; i < output.size(); ++i) {
            const auto value = table.powers[j][i];
            output[i] = f.add(output[i], poly[j] == 2 ? f.negate(value) : value);
        }
    }
}

#ifdef CENSUS_AVX2
__attribute__((target("avx2")))
__m256i add3(__m256i a, __m256i b) {
    const auto mask = _mm256_set1_epi32(65535);
    const auto al = _mm256_and_si256(a, mask), ah = _mm256_srli_epi32(a, 16);
    const auto bl = _mm256_and_si256(b, mask), bh = _mm256_srli_epi32(b, 16);
    const auto an = _mm256_or_si256(al, ah), bn = _mm256_or_si256(bl, bh);
    const auto lo = _mm256_or_si256(_mm256_or_si256(
        _mm256_andnot_si256(bn, al), _mm256_andnot_si256(an, bl)), _mm256_and_si256(ah, bh));
    const auto hi = _mm256_or_si256(_mm256_or_si256(
        _mm256_andnot_si256(bn, ah), _mm256_andnot_si256(an, bh)), _mm256_and_si256(al, bl));
    return _mm256_or_si256(lo, _mm256_slli_epi32(hi, 16));
}

__attribute__((target("avx2")))
void evaluate_avx2(const CountingTable& table, const std::vector<int>& poly,
                   std::vector<std::uint32_t>& output) {
    output.assign(table.powers[0].size(), 0);
    const auto& f = table.traces.field;
    for (std::size_t j = 0; j < poly.size(); ++j) {
        if (!poly[j]) continue;
        std::size_t i = 0;
        for (; i + 8 <= output.size(); i += 8) {
            auto v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(table.powers[j].data()+i));
            if (poly[j] == 2) v = _mm256_or_si256(_mm256_srli_epi32(v, 16), _mm256_slli_epi32(v, 16));
            const auto a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output.data()+i));
            const auto sum = f.p == 2 ? _mm256_xor_si256(a, v) : add3(a, v);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(output.data()+i), sum);
        }
        for (; i < output.size(); ++i) {
            const auto value = table.powers[j][i];
            output[i] = f.add(output[i], poly[j] == 2 ? f.negate(value) : value);
        }
    }
}

__attribute__((target("avx2")))
bool sum_ordinary3(const CountingTable& table, const std::vector<std::uint32_t>& c,
                   const std::vector<std::uint32_t>& delta, std::int64_t& sum) {
    const auto& f = table.traces.field;
    const auto zero = _mm256_setzero_si256(), mask = _mm256_set1_epi32(f.mask);
    const auto six = _mm256_set1_epi32(6), sixm = _mm256_set1_epi32(6*f.order);
    sum = 0;
    std::size_t i = 0;
    for (; i + 8 <= c.size(); i += 8) {
        const auto a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(c.data()+i));
        const auto b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(delta.data()+i));
        if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmpeq_epi32(a, zero), _mm256_cmpeq_epi32(b, zero)))) return false;
        const auto index = [&](const __m256i v) __attribute__((target("avx2"))) {
            return _mm256_add_epi32(_mm256_i32gather_epi32(
                reinterpret_cast<const int*>(f.index_mask.data()), _mm256_and_si256(v, mask), 4),
                _mm256_slli_epi32(_mm256_i32gather_epi32(
                reinterpret_cast<const int*>(f.index_mask.data()), _mm256_srli_epi32(v, 16), 4), 1));
        };
        const auto la = _mm256_i32gather_epi32(f.logarithms.data(), index(a), 4);
        const auto lb = _mm256_i32gather_epi32(f.logarithms.data(), index(b), 4);
        const auto e = _mm256_add_epi32(_mm256_sub_epi32(lb, _mm256_mullo_epi32(la, six)), sixm);
        alignas(32) int exponents[8];
        _mm256_store_si256(reinterpret_cast<__m256i*>(exponents), e);
        for (int& exponent : exponents) exponent %= f.order;
        auto tr = _mm256_i32gather_epi32(table.traces.ordinary.data(),
            _mm256_load_si256(reinterpret_cast<const __m256i*>(exponents)), 4);
        const auto sign = _mm256_sub_epi32(zero, _mm256_and_si256(la, _mm256_set1_epi32(1)));
        tr = _mm256_sub_epi32(_mm256_xor_si256(tr, sign), sign);
        alignas(32) int values[8];
        _mm256_store_si256(reinterpret_cast<__m256i*>(values), tr);
        for (const int value : values) sum += value;
    }
    for (; i < c.size(); ++i) {
        if (!c[i] || !delta[i]) return false;
        const int a = f.logarithms[f.index(c[i])], b = f.logarithms[f.index(delta[i])];
        sum += (a & 1 ? -1 : 1)*table.traces.ordinary[(b-6*a+6*f.order) % f.order];
    }
    return true;
}
#endif

std::uint32_t power(const TraceField& f, std::uint32_t a, int n) {
    return a ? f.exp(std::int64_t(f.logarithms[f.index(a)])*n) : 0;
}

std::uint32_t divide(const TraceField& f, std::uint32_t a, std::uint32_t b) {
    return a ? f.exp(f.logarithms[f.index(a)]-f.logarithms[f.index(b)]) : 0;
}

int supersingular2(const TraceField& f, std::uint32_t a2, std::uint32_t a3,
                   std::uint32_t a4, std::uint32_t a6) {
    const auto inverse = f.exp(-2LL*f.logarithms[f.index(a3)]);
    const auto u = inverse, v = f.multiply(a2, inverse);
    const auto w = f.multiply(a4, inverse), z = f.multiply(a6, inverse);
    unsigned linear = 0, active = (1u << f.degree)-1;
    unsigned cross[12]{};
    int constant = f.trace(z), factor = 1;
    for (int i = 0; i < f.degree; ++i) {
        const auto x = 1u << i;
        linear |= unsigned(f.trace(f.add(f.add(f.multiply(u, power(f, x, 3)),
            f.multiply(v, power(f, x, 2))), f.multiply(w, x)))) << i;
        for (int j = 0; j < i; ++j) {
            const auto y = 1u << j;
            const auto b = f.trace(f.multiply(u, f.add(f.multiply(power(f, x, 2), y),
                f.multiply(x, power(f, y, 2)))));
            if (b) { cross[i] |= 1u << j; cross[j] |= 1u << i; }
        }
    }
    // Eliminate a hyperbolic pair in Tr((x^3+a2*x^2+a4*x+a6)/a3^2).
    while (true) {
        int i = -1, j = -1;
        for (int k = 0; k < f.degree; ++k) {
            if ((active >> k & 1) && (cross[k] & active)) {
                i = k;
                j = __builtin_ctz(cross[k] & active);
                break;
            }
        }
        if (i < 0) break;
        active &= ~((1u << i) | (1u << j));
        const unsigned a = cross[i] & active, b = cross[j] & active;
        const int ai = linear >> i & 1, bj = linear >> j & 1;
        constant ^= ai & bj;
        linear ^= (ai ? b : 0) ^ (bj ? a : 0) ^ (a & b);
        for (int k = 0; k < f.degree; ++k) {
            if (active >> k & 1) {
                const unsigned row = ((a >> k & 1) ? b : 0) ^ ((b >> k & 1) ? a : 0);
                cross[k] ^= row & ~(1u << k);
            }
        }
        factor *= 2;
    }
    if (linear & active) return 0;
    return -(constant ? -1 : 1)*factor*(1 << __builtin_popcount(active));
}

template<class Values>
int good_trace(const CountingTable& table, std::uint32_t c, std::uint32_t delta, Values value) {
    const auto& f = table.traces.field;
    if (f.p == 3) {
        if (c) {
            const int log = f.logarithms[f.index(c)];
            const int e = (f.logarithms[f.index(delta)]-6*log+6*f.order) % f.order;
            return (log & 1 ? -1 : 1)*table.traces.ordinary[e];
        }
        const auto a = value(3), b = value(4);
        const int e = f.logarithms[f.index(f.negate(a))];
        if (e & 1) return 0;
        const auto r = f.exp(e/2);
        const int sign = f.logarithms[f.index(r)] & 1 ? -1 : 1;
        return -3*sign*table.gauss[f.trace(divide(f, b, power(f, r, 3)))];
    }
    const auto a2 = value(1), a3 = value(2);
    if (!c) return supersingular2(f, a2, a3, value(3), value(4));
    const auto twist = f.add(divide(f, a2, power(f, c, 2)), divide(f, a3, power(f, c, 3)));
    const int e = (f.logarithms[f.index(delta)]-12*f.logarithms[f.index(c)]+12*f.order) % f.order;
    return (f.trace(twist) ? -1 : 1)*table.traces.ordinary[e];
}

void local_series(std::vector<std::int64_t>& sums, int trace, bool good, int norm, int degree) {
    std::int64_t previous = 2, current = trace;
    for (std::size_t n = degree; n <= sums.size(); n += degree) {
        sums[n-1] += degree*current;
        const auto next = trace*current-(good ? norm*previous : 0);
        previous = current;
        current = next;
    }
}
}

struct LPolynomialCounter::State {
    bool avx2 = false;
    std::map<std::pair<int, int>, CountingTable> tables;

    CountingTable& table(int p, int d) {
        return tables.try_emplace({p, d}, p, d).first->second;
    }
};

LPolynomialCounter::LPolynomialCounter(std::string_view backend) : state(std::make_unique<State>()) {
    if (backend != "auto" && backend != "scalar" && backend != "avx2") {
        throw std::invalid_argument("backend must be auto, scalar or avx2");
    }
#ifdef CENSUS_AVX2
    state->avx2 = backend != "scalar" && __builtin_cpu_supports("avx2");
#endif
    if (backend == "avx2" && !state->avx2) {
        throw std::invalid_argument("AVX2 backend is unavailable");
    }
}

LPolynomialCounter::~LPolynomialCounter() = default;

LPolynomial LPolynomialCounter::compute(const Model& model, const LocalModels& local,
                                        const SignData& sign) {
    NTL::zz_pPush context(model.q);
    const auto a = polynomials(model.affine);
    const auto c = model.q == 2 ? a[0] : a[1]+a[0]*a[0];
    const auto discriminant = WeierstrassInvariants(a).discriminant();
    const auto delta = model.q == 2 ? discriminant : -discriminant;
    const auto coefficients = [](const NTL::zz_pX& poly) {
        std::vector<int> result(NTL::deg(poly)+1);
        for (long j = 0; j <= NTL::deg(poly); ++j) result[j] = NTL::rep(NTL::coeff(poly, j));
        return result;
    };
    const auto primary_poly = coefficients(c), delta_poly = coefficients(delta);
    LPolynomial result;
    result.trace_sums.resize(sign.degree/2);
    for (const auto& fibre : local.fibres) {
        local_series(result.trace_sums, fibre.trace, false, 0, fibre.degree);
    }
    for (int d = 1; d <= sign.degree/2; ++d) {
        const auto& table = state->table(model.q, d);
        const auto& f = table.traces.field;
        std::vector<std::uint32_t> primary, values;
#ifdef CENSUS_AVX2
        const auto evaluate = state->avx2 ? evaluate_avx2 : evaluate_scalar;
#else
        const auto evaluate = evaluate_scalar;
#endif
        evaluate(table, primary_poly, primary);
        evaluate(table, delta_poly, values);
        bool summed = false;
#ifdef CENSUS_AVX2
        std::int64_t sum;
        if (state->avx2 && model.q == 3 && d > sign.degree/4
            && sum_ordinary3(table, primary, values, sum)) {
            result.trace_sums[d-1] += d*sum;
            summed = true;
        }
#endif
        if (!summed) {
            for (std::size_t i = 0; i < primary.size(); ++i) {
                if (!values[i]) continue;
                const auto value = [&](int j) {
                    std::uint32_t v = 0;
                    for (std::size_t k = 0; k < model.affine[j].size(); ++k) {
                        const int digit = model.affine[j][k];
                        if (digit) v = f.add(v, digit == 2
                            ? f.negate(table.powers[k][i]) : table.powers[k][i]);
                    }
                    return v;
                };
                local_series(result.trace_sums, good_trace(table, primary[i], values[i], value),
                    true, f.size, d);
            }
        }
        if (d == 1 && std::none_of(local.fibres.begin(), local.fibres.end(),
                                   [](const auto& fibre) { return fibre.place == "inf"; })) {
            const auto value = [&](int j) {
                const auto digit = model.infinity[j][0];
                return digit == 2 ? f.negate(1) : std::uint32_t(digit);
            };
            const auto c_infinity = model.q == 2 ? value(0) : f.add(value(1), power(f, value(0), 2));
            const int digit = NTL::rep(NTL::coeff(delta, 24));
            const auto delta_infinity = digit == 2 ? f.negate(1) : std::uint32_t(digit);
            local_series(result.trace_sums, good_trace(table, c_infinity, delta_infinity, value),
                true, model.q, 1);
        }
    }
    result.coefficients.resize(sign.degree+1);
    result.coefficients[0] = 1;
    for (int n = 1; n <= sign.degree/2; ++n) {
        boost::multiprecision::cpp_int numerator = 0;
        for (int j = 1; j <= n; ++j) numerator += result.trace_sums[j-1]*result.coefficients[n-j];
        if (numerator % n != 0) throw std::runtime_error("nonintegral Euler coefficient");
        result.coefficients[n] = numerator/n;
    }
    for (int i = 0; i <= sign.degree/2; ++i) {
        boost::multiprecision::cpp_int factor = sign.sign;
        for (int j = 0; j < sign.degree-2*i; ++j) factor *= model.q;
        const boost::multiprecision::cpp_int mirrored = factor*result.coefficients[i];
        if (sign.degree-i == i && mirrored != result.coefficients[i]) {
            throw std::runtime_error("functional equation middle coefficient is nonzero");
        }
        result.coefficients[sign.degree-i] = mirrored;
    }
    return result;
}

}
