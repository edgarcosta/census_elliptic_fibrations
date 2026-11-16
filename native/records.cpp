#include "records.hpp"

#include <sstream>
#include <stdexcept>
#include <zstd.h>

namespace census {
namespace {
std::string quote(std::string_view cell) {
    std::string result;
    const bool quoted = cell.find_first_of(",\"\r\n") != cell.npos;
    if (quoted) result += '"';
    for (const unsigned char c : cell) {
        if (c == 0 || c > 127) throw std::invalid_argument("record cells must be ASCII without NUL");
        if (c == '"') result += '"';
        result += c;
    }
    if (quoted) result += '"';
    return result;
}

template<class Values>
std::string wire_list(const Values& values) {
    if (values.empty()) return "[]";
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ' ';
        out << values[i];
    }
    return out.str();
}

template<class Values>
std::string optional_list(const std::optional<Values>& values) {
    return values ? wire_list(*values) : "";
}

void check_zstd(std::size_t result) {
    if (ZSTD_isError(result)) throw std::runtime_error(ZSTD_getErrorName(result));
}
}

InvariantRow invariant_row(const Model& model, LPolynomialCounter& counter) {
    InvariantRow row{};
    row.q = model.q;
    row.key = model.key;
    row.kind = model.kind;
    row.minimality = model.minimality;
    row.status = model.status;
    if (model.status == "unsupported") {
        row.refinement_status = "unsupported";
        return row;
    }
    try {
        const auto local = local_models(model);
        row.fibres = local.fibres;
        const auto sign = functional_equation_sign(local, model.q);
        row.degree = sign.degree;
        row.sign = sign.sign;
        const auto torsion = complete_torsion(model, local);
        row.torsion = torsion.factors;
        const auto lfun = counter.compute(model, local, sign);
        row.lower_l = std::vector<boost::multiprecision::cpp_int>(
            lfun.coefficients.begin()+1, lfun.coefficients.begin()+1+sign.degree/2);
        const auto derived = derived_invariants(model, local, torsion, lfun);
        row.sha_reg = derived.sha_reg;
        const auto candidates = admissible_candidates(derived, local.fibres);
        row.admissible_sha = candidates.sha;
        const bool singleton = candidates.sha.size() == 1;
        row.split_method = singleton ? "constraints" : "unresolved";
        row.refinement_status = singleton ? "not_needed" : model.q == 2 ? "unsupported" : "";
        row.status = candidates.sha.empty() ? "error" : "ok";
        if (candidates.sha.empty()) row.message = "no admissible Sha order";
    } catch (const std::exception& error) {
        row.status = "error";
        row.refinement_status = "error";
        row.message = error.what();
    }
    return row;
}

std::string invariant_csv_header() {
    return "schema_version,q,key,kind,minimality,status,degree,sign,lower_L,fibres,torsion,"
           "sha_reg,admissible_sha,split_method,d2,sha2_dimensions,refinement_status,message";
}

std::string invariant_csv(const InvariantRow& row) {
    std::string fibres;
    if (row.fibres) {
        if (row.fibres->empty()) fibres = "[]";
        else {
            std::ostringstream out;
            for (std::size_t i = 0; i < row.fibres->size(); ++i) {
                if (i) out << ';';
                const auto& f = (*row.fibres)[i];
                out << f.place << ':' << f.degree << ':' << f.type << ':' << f.index
                    << ':' << f.discriminant_valuation << ':' << f.conductor << ':'
                    << f.tamagawa << ':' << f.trace;
            }
            fibres = out.str();
        }
    }
    std::string rational;
    if (row.sha_reg) {
        std::ostringstream out;
        out << numerator(*row.sha_reg) << '/' << denominator(*row.sha_reg);
        rational = out.str();
    }
    const std::array<std::string, 18> cells = {
        "1", std::to_string(row.q), std::to_string(row.key), row.kind, row.minimality, row.status,
        row.degree ? std::to_string(*row.degree) : "", row.sign ? std::to_string(*row.sign) : "",
        optional_list(row.lower_l), fibres, optional_list(row.torsion), rational,
        optional_list(row.admissible_sha), row.split_method, "", optional_list(row.sha2_dimensions),
        row.refinement_status, row.message};
    std::string text;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        if (i) text += ',';
        text += quote(cells[i]);
    }
    return text;
}

struct CompressedRecords::State {
    std::ostream& output;
    ZSTD_CCtx* context = ZSTD_createCCtx();
    std::vector<char> buffer = std::vector<char>(ZSTD_CStreamOutSize());
    bool finished = false;

    explicit State(std::ostream& out) : output(out) {
        if (!context) throw std::runtime_error("could not create zstd compressor");
    }
    ~State() { ZSTD_freeCCtx(context); }

    void compress(std::string_view text, ZSTD_EndDirective directive) {
        ZSTD_inBuffer input{text.data(), text.size(), 0};
        std::size_t remaining;
        do {
            ZSTD_outBuffer out{buffer.data(), buffer.size(), 0};
            remaining = ZSTD_compressStream2(context, &out, &input, directive);
            check_zstd(remaining);
            output.write(buffer.data(), out.pos);
            if (!output) throw std::runtime_error("could not write standard output");
        } while (input.pos < input.size || (directive == ZSTD_e_end && remaining));
    }
};

CompressedRecords::CompressedRecords(std::ostream& output) : state(std::make_unique<State>(output)) {
    check_zstd(ZSTD_CCtx_setParameter(state->context, ZSTD_c_compressionLevel, 9));
    check_zstd(ZSTD_CCtx_setParameter(state->context, ZSTD_c_checksumFlag, 1));
}
CompressedRecords::~CompressedRecords() = default;
void CompressedRecords::write(std::string_view text) { state->compress(text, ZSTD_e_continue); }
void CompressedRecords::finish() {
    if (state->finished) return;
    state->compress({}, ZSTD_e_end);
    state->output.flush();
    if (!state->output) throw std::runtime_error("could not write standard output");
    state->finished = true;
}
}
