import csv
from fractions import Fraction
import io
import math
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from _support import ROOT, invoke


FIELDS = ("admissible_sha", "admissible_regulators")


def values(text, convert=int):
    return [] if text == "[]" else [convert(v) for v in text.split()]


def check_candidates(row, want):
    for field in FIELDS:
        assert row.get(field) == want[field], (field, row.get(field), want[field], want)
    sha = values(row["admissible_sha"])
    regulators = values(row["admissible_regulators"], Fraction)
    assert sha == sorted(set(sha))
    assert all(s > 0 and math.isqrt(s)**2 == s for s in sha)
    assert len(sha) == len(regulators)
    assert all(s*r == Fraction(want["sha_reg"]) for s, r in zip(sha, regulators))


def test_equations(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert len(expected) == 16
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    result = invoke(text)
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected)
    for row, want in zip(actual, expected):
        assert (row["q"], row["key"]) == (want["q"], want["key"])
        for field in FIELDS:
            assert row[field] == want[field], (field, row, want)
    assert any(r["rank"] == "0" for r in actual)
    assert any(r["admissible_sha"] == "1 4 9" for r in actual)
    assert any(r["admissible_regulators"] == "14 7/2" for r in actual)
    assert any(r["admissible_regulators"] == "3/4" for r in actual)


# Constant extensions and large supplied integers have no prime-field census key.
def test_supplied(reference):
    result = subprocess.run([ROOT / "build/tests/native/candidates-exact"],
                            input=Path(reference).read_text(), text=True,
                            capture_output=True)
    assert result.returncode == 0, result.stderr
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(expected) == 192
    assert len(actual) == len(expected)
    for row, want in zip(actual, expected):
        check_candidates(row, want)
    assert any(max(values(r["admissible_sha"]), default=0) > 2**53 for r in actual)
    assert any(r["admissible_sha"] == "1 9 81 729" for r in actual)
    assert any(r["admissible_sha"] == "[]" for r in actual)


if __name__ == "__main__":
    test_equations(sys.argv[1])
    test_supplied(sys.argv[1] + ".supplied")
