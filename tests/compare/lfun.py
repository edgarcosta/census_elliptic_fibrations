import csv
import io
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from _support import invoke


MIDDLE_ZERO_SIGNS = {(2, 4884480): 1, (2, 11821184): -1,
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
        q, degree, sign = int(got["q"]), int(got["degree"]), int(got["sign"])
        coefficients = integers(got["L_polynomial"])
        if (q, int(got["key"])) in MIDDLE_ZERO_SIGNS:
            assert sign == MIDDLE_ZERO_SIGNS[q, int(got["key"])], context
            assert all(c == 0 for c in coefficients[degree // 2:13]), context


# A wrong Euler factor, chart, sign, or Newton coefficient changes this result.
def test_magma_comparison(reference):
    text = "".join(f"{r['q']} {r['key']}\n" for r in reference)
    assert MIDDLE_ZERO_SIGNS.keys() <= {(int(r["q"]), int(r["key"])) for r in reference}
    automatic = invoke(text)
    assert automatic.returncode == 0, automatic.stderr
    compare(reference, automatic.stdout)
    scalar = invoke(text, "--backend", "scalar")
    assert scalar.returncode == 0, scalar.stderr
    assert scalar.stdout == automatic.stdout
    avx2 = invoke(text, "--backend", "avx2")
    if avx2.returncode:
        assert "AVX2 backend is unavailable" in avx2.stderr, avx2.stderr
    else:
        assert avx2.stdout == scalar.stdout



if __name__ == "__main__":
    with open(sys.argv[1], newline="") as source:
        reference = list(csv.DictReader(source))
    assert reference
    test_magma_comparison(reference)
