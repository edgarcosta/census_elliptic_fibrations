#include "records.hpp"
#include "format.hpp"

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

void check_zstd(std::size_t result) {
    if (ZSTD_isError(result)) throw std::runtime_error(ZSTD_getErrorName(result));
}
}

std::string invariant_csv_header() {
    return "schema_version,q,key,kind,minimality,status,degree,sign,lower_L,fibres,torsion,"
           "sha_reg,admissible_sha,split_method,d2,sha2_dimensions,refinement_status,message";
}

std::string invariant_csv(const Model& model, const InvariantData& data,
    const std::optional<std::string>& error) {
    std::array<std::string, 18> cells = {};
    cells[0] = "1";
    cells[1] = std::to_string(model.q);
    cells[2] = std::to_string(model.key);
    cells[3] = model.kind;
    cells[4] = model.minimality;
    cells[5] = model.status;
    if (model.status == "unsupported") cells[16] = "unsupported";
    if (data.local) {
        std::ostringstream out;
        write_fibres(out, data.local->fibres);
        cells[9] = out.str();
    }
    if (data.sign) {
        cells[6] = std::to_string(data.sign->degree);
        cells[7] = std::to_string(data.sign->sign);
    }
    if (data.torsion) {
        std::ostringstream out;
        write_list(out, data.torsion->factors);
        cells[10] = out.str();
    }
    if (data.lfun) {
        std::ostringstream out;
        write_list(out, data.lfun->coefficients, 1, 1+data.sign->degree/2);
        cells[8] = out.str();
    }
    if (data.derived) {
        std::ostringstream out;
        out << numerator(data.derived->sha_reg) << '/' << denominator(data.derived->sha_reg);
        cells[11] = out.str();
    }
    if (data.candidates) {
        std::ostringstream out;
        write_list(out, data.candidates->sha);
        cells[12] = out.str();
        cells[13] = data.candidates->split_method();
        cells[16] = cells[13] == "constraints" ? "not_needed" : model.q == 2 ? "unsupported" : "";
        cells[5] = data.candidates->sha.empty() ? "error" : "ok";
        if (data.candidates->sha.empty()) cells[17] = "no admissible Sha order";
    }
    if (error) {
        cells[5] = "error";
        cells[16] = "error";
        cells[17] = *error;
    }
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
