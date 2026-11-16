import csv
from fractions import Fraction
import io
import math
import os
from pathlib import Path
import subprocess
import sys
from _support import fake_magma


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("CENSUS_NATIVE", ROOT / "artifacts/native/census-native"))
FIELDS = ("admissible_sha", "admissible_regulators", "split_method",
          "regulator", "sha", "brauer_order", "ns_discriminant")


def values(text, convert=int):
    return [] if text == "[]" else [convert(v) for v in text.split()]


def check_candidates(row, want):
    for field in FIELDS[:3]:
        assert row.get(field) == want[field], (field, row.get(field), want[field], want)
    sha = values(row["admissible_sha"])
    regulators = values(row["admissible_regulators"], Fraction)
    assert sha == sorted(set(sha))
    assert all(s > 0 and math.isqrt(s)**2 == s for s in sha)
    assert len(sha) == len(regulators)
    assert all(s*r == Fraction(want["sha_reg"]) for s, r in zip(sha, regulators))
    if len(sha) == 1:
        assert row["split_method"] == "constraints"
        assert Fraction(row["regulator"]) == regulators[0]
        assert int(row["sha"]) == int(row["brauer_order"]) == sha[0]
        assert int(row["ns_discriminant"]) * sha[0] == int(want["artin_tate_product"])
    else:
        assert row["split_method"] == "unresolved"
        assert all(row[f] == "" for f in FIELDS[3:])


def test_equations(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert len(expected) == 16
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    with fake_magma() as (env, marker):
        result = subprocess.run([BINARY, "invariants", "--diagnostic"], input=text, text=True,
                                capture_output=True, env=env, cwd=ROOT)
        assert not marker.exists(), "native command started Magma"
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected)
    for row, want in zip(actual, expected):
        assert (row["q"], row["key"]) == (want["q"], want["key"])
        check_candidates(row, want)
    assert any(r["rank"] == "0" for r in actual)
    assert any(r["admissible_sha"] == "1 4 9" for r in actual)
    assert any(r["admissible_regulators"] == "14 7/2" for r in actual)
    assert any(r["regulator"] == "3/4" for r in actual)
    print(f"CANDIDATES_NATIVE_OK: {len(actual)} records, complete lists and proved splits")


# Constant extensions and large supplied integers have no prime-field census key.
def test_supplied(reference):
    result = subprocess.run([ROOT / "artifacts/native/candidates-exact"],
                            input=Path(reference).read_text(), text=True,
                            capture_output=True, check=True)
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(expected) == 192
    assert len(actual) == len(expected)
    for row, want in zip(actual, expected):
        want["split_method"] = "constraints" if len(values(want["admissible_sha"])) == 1 else "unresolved"
        check_candidates(row, want)
    assert any(max(values(r["admissible_sha"]), default=0) > 2**53 for r in actual)
    assert any(r["admissible_sha"] == "1 9 81 729" for r in actual)
    assert any(r["admissible_sha"] == "[]" for r in actual)
    print(f"CANDIDATES_EXACT_OK: {len(actual)} supplied controls equal Magma")


def test_unsupported():
    result = subprocess.run([BINARY, "invariants", "--diagnostic"], input="2 0\n3 0\n2 1\n3 1\n",
                            text=True, capture_output=True, check=True)
    for row in csv.DictReader(io.StringIO(result.stdout)):
        assert row["status"] == "unsupported"
        assert all(row[f] == "" for f in FIELDS)


if __name__ == "__main__":
    test_equations(sys.argv[1])
    test_supplied(sys.argv[1] + ".supplied")
    test_unsupported()
