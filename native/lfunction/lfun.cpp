#include "lfun.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

namespace census {

SignData functional_equation_sign(const LocalModels& local, unsigned q) {
    auto summarize = [q](const std::vector<LocalFibre>& fibres) {
        std::array<int, 3> values{};
        for (const auto& fibre : fibres) {
            values[0] += fibre.degree*fibre.discriminant_valuation;
            values[1] += fibre.degree*fibre.conductor;
            // Each closed place contributes once to C, without a degree exponent.
            for (int c = fibre.tamagawa; c % q == 0; c /= q) ++values[2];
        }
        if (values[0] % 12) throw std::runtime_error("nonintegral elliptic surface height");
        return values;
    };
    const auto original = summarize(local.fibres);
    const auto twist = summarize(local.twist_fibres);
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
    // M > 2Q makes centered recovery unique in [-Q,Q]; see docs/native/tables.md.
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


}

TraceTable::TraceTable(int p, int degree)
    : field(p, degree), nodal(0) {
    if (p == 3) ordinary3(*this);
    else ordinary2(*this);
    if (nodal != 1) throw std::runtime_error("nodal trace is not one");
}

}
