import csv
import io
import os
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from _support import ROOT, fake_magma, invoke


LFUN_FIELDS = {"degree", "sign", "lower_L", "L_polynomial", "trace_sums"}
FORMER_OPEN = {(2, 4884480): 1, (2, 11821184): -1,
               (3, 372699133912): -1, (3, 320009441124): 1}


def integers(text):
    return [] if text == "[]" else list(map(int, text.split()))


def check_identities(row):
    context = (row["q"], row["key"])
    q, degree, sign = int(row["q"]), int(row["degree"]), int(row["sign"])
    assert sign in (-1, 1), context
    coefficients = integers(row["L_polynomial"])
    assert len(coefficients) == degree + 1 and coefficients[0] == 1, context
    assert integers(row["lower_L"]) == coefficients[1:1 + degree // 2], context
    for i in range(degree // 2 + 1):
        assert coefficients[degree - i] == sign * q ** (degree - 2*i) * coefficients[i], context


def compare(reference, output):
    actual = list(csv.DictReader(io.StringIO(output)))
    assert len(actual) == len(reference)
    for got, want in zip(actual, reference):
        context = (want["q"], want["key"])
        for field, value in want.items():
            assert got[field] == value, (context, field, got.get(field), value)
        check_identities(got)
        q, degree, sign = int(got["q"]), int(got["degree"]), int(got["sign"])
        coefficients = integers(got["L_polynomial"])
        if (q, int(got["key"])) in FORMER_OPEN:
            assert sign == FORMER_OPEN[q, int(got["key"])], context
            assert all(c == 0 for c in coefficients[degree // 2:13]), context
    return actual


# A wrong Euler factor, chart, sign, or Newton coefficient changes this result.
def test_magma_comparison(reference):
    text = "".join(f"{r['q']} {r['key']}\n" for r in reference)
    assert FORMER_OPEN.keys() <= {(int(r["q"]), int(r["key"])) for r in reference}
    with fake_magma() as (env, marker):
        automatic = invoke(text, env=env)
        assert automatic.returncode == 0, automatic.stderr
        assert not marker.exists(), "native command started Magma"
    rows = compare(reference, automatic.stdout)
    scalar = invoke(text, "--backend", "scalar")
    assert scalar.returncode == 0, scalar.stderr
    assert scalar.stdout == automatic.stdout
    avx2 = invoke(text, "--backend", "avx2")
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists() and "avx2" in cpuinfo.read_text().split():
        assert avx2.returncode == 0, avx2.stderr
    if avx2.returncode:
        assert "AVX2 backend is unavailable" in avx2.stderr, avx2.stderr
    else:
        assert avx2.stdout == scalar.stdout
    if "CENSUS_NATIVE" not in os.environ:
        scalar_binary = ROOT / "artifacts/native/census-native-scalar"
        result = subprocess.run([scalar_binary, "invariants", "--diagnostic"], input=text,
                                capture_output=True, text=True, cwd=ROOT)
        assert result.returncode == 0, result.stderr
        assert result.stdout == scalar.stdout
        unavailable = subprocess.run([scalar_binary, "invariants", "--diagnostic", "--backend", "avx2"],
                                     input=text, capture_output=True, text=True, cwd=ROOT)
        assert unavailable.returncode != 0 and "AVX2 backend is unavailable" in unavailable.stderr
    return rows


# Coefficients determine the key and sign without trusting caller-supplied data.
def test_coefficient_input(reference, rows):
    def code(q, coefficients):
        return sum(c * q**i for i, c in enumerate(integers(coefficients)))

    text = "".join(" ".join([row["q"], *[str(code(int(row["q"]), row[a]))
                   for a in ("a1", "a2", "a3", "a4", "a6")]]) + "\n" for row in rows)
    result = invoke(text, "--coefficients", "--backend", "scalar")
    assert result.returncode == 0, result.stderr
    compare(reference, result.stdout)
    row = text.splitlines()[0]
    for text, diagnostic in (("2 0 1 1 0 0\n", "census normal form"),
                             ("3 1 0 0 0 0\n", "census normal form"),
                             ("3 0 0 1 0 0\n", "census normal form"),
                             ("2 8 0 1 0 0\n", "K3 degree bound")):
        invalid = invoke(text, "--coefficients")
        assert invalid.returncode != 0 and diagnostic in invalid.stderr, invalid
        assert not list(csv.DictReader(io.StringIO(invalid.stdout))), invalid.stdout
    supplied = invoke(row + " 1\n", "--coefficients")
    assert supplied.returncode != 0 and "expected q a1 a2 a3 a4 a6" in supplied.stderr


def test_unsupported_cells():
    result = invoke("2 0\n3 0\n2 1\n3 1\n")
    assert result.returncode == 0, result.stderr
    for row in csv.DictReader(io.StringIO(result.stdout)):
        assert row["status"] == "unsupported", row
        assert all(row[field] == "" for field in LFUN_FIELDS), row


if __name__ == "__main__":
    with open(sys.argv[1], newline="") as source:
        reference = list(csv.DictReader(source))
    assert reference
    rows = test_magma_comparison(reference)
    test_coefficient_input(reference, rows)
    test_unsupported_cells()
