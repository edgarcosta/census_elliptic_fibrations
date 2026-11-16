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
FIELDS = ("rank", "geometric_rank", "picard_number", "geometric_picard_number",
          "P2", "zeta_sign", "point_counts", "brauer_height", "artin_tate_product",
          "ns_splitting_degree", "j_degree", "j_inseparability_exponent",
          "is_isotrivial", "sha_reg")


def invoke(text, **kwargs):
    return subprocess.run([BINARY, "invariants"], input=text, capture_output=True,
                          text=True, cwd=ROOT, **kwargs)


def integers(text):
    return [] if text == "[]" else list(map(int, text.split()))


# Closed-place Tamagawa factors occur once; cycles already include residue degree.
def check_identities(row):
    q = int(row["q"])
    fibres = [f.split(":") for f in row["fibres"].split(";")]
    cycles = [n for f in fibres for n in integers(f[8])]
    rank = int(row["rank"])
    assert int(row["picard_number"]) == 2 + rank + len(cycles), row
    assert int(row["geometric_picard_number"]) == 2 + int(row["geometric_rank"]) + sum(cycles), row
    p2 = integers(row["P2"])
    assert len(p2) == 23 and p2[-1] == int(row["zeta_sign"]) * q**22, row
    traces = [0]
    for n in range(1, 6):
        traces.append(-n*p2[n] - sum(p2[j]*traces[n-j] for j in range(1, n)))
    assert integers(row["point_counts"]) == [1+q**(2*n)+traces[n] for n in range(1, 6)], row
    torsion = math.prod(integers(row["torsion"]))
    tamagawa = math.prod(int(f[6]) for f in fibres)
    assert Fraction(row["artin_tate_product"]) == Fraction(row["sha_reg"]) * tamagawa * math.prod(cycles) / torsion**2, row
    assert (row["is_isotrivial"] == "true") == (int(row["j_degree"]) == 0), row
    assert int(row["ns_splitting_degree"]) % math.lcm(1, *cycles) == 0, row


def test_magma_comparison(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert len(expected) == 28, len(expected)
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    with fake_magma() as (env, marker):
        result = invoke(text, env=env)
        assert not marker.exists(), "native command started Magma"
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        for field in FIELDS:
            assert got.get(field) == want[field], (want["q"], want["key"], field, got.get(field), want[field])
        check_identities(got)
    assert {r["brauer_height"] for r in actual} >= {"1", "3", "4", "8", "infinity"}
    assert any(r["rank"] == "0" and r["sha_reg"] == "1" and r["torsion"] == "6" for r in actual)
    assert any(r["sha_reg"] == "3/4" for r in actual)
    assert {r["j_inseparability_exponent"] for r in actual} >= {"0", "1"}
    assert any(int(f.split(":")[1]) > 1 and int(f.split(":")[6]) > 1
               for r in actual for f in r["fibres"].split(";"))
    print(f"DERIVED_NATIVE_OK: {len(actual)} records, {len(FIELDS)} fields and exact identities")


def test_unsupported():
    result = invoke("2 0\n3 0\n2 1\n3 1\n")
    assert result.returncode == 0, result.stderr
    for row in csv.DictReader(io.StringIO(result.stdout)):
        assert row["status"] == "unsupported", row
        assert all(row[field] == "" for field in FIELDS), row


# Prime-field keys cannot exercise q=4096; the supplied-data API uses the same engine.
def test_constant_field_controls(reference):
    exact = subprocess.run([ROOT / "artifacts/native/derived-exact"], capture_output=True, text=True, check=True)
    assert exact.stdout == Path(reference + ".exact").read_text()
    large = dict(zip(FIELDS, exact.stdout.splitlines()[0].split(",")))
    assert integers(large["P2"])[-1] == 2**264
    assert integers(large["point_counts"])[-1] == 1329227995784915898268080161630978049
    print("DERIVED_EXACT_OK: GF(4) and GF(4096), all fields equal Magma")


if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
    test_unsupported()
    test_constant_field_controls(sys.argv[1])
