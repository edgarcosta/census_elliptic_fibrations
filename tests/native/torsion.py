import csv
import io
from pathlib import Path
import subprocess
import sys
from _support import fake_magma


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"


def invoke(text, **kwargs):
    return subprocess.run([BINARY, "invariants"], input=text, capture_output=True,
                          text=True, cwd=ROOT, **kwargs)


def test_magma_comparison(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert len(expected) == 26, len(expected)
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    with fake_magma() as (env, marker):
        result = invoke(text, env=env)
        assert not marker.exists(), "native command started Magma"
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        context = (want["q"], want["key"])
        assert got["torsion"] == want["torsion"], context
        points = got["torsion_points"].split(";")
        orders = list(map(int, got["torsion_orders"].split()))
        wanted = dict(zip(want["torsion_points"].split(";"),
                          map(int, want["torsion_orders"].split())))
        assert len(points) == len(set(points)) == len(orders), context
        assert dict(zip(points, orders)) == wanted, context
        assert points[0] == "O" and orders[0] == 1, context
    assert {r["torsion"] for r in actual} >= {"[]", "2", "3", "6", "2 4"}
    assert any("1:4:1:8" in r["torsion_points"].split(";") for r in actual), "degree-one denominator"
    assert any("64:1:512:1" in r["torsion_points"].split(";") for r in actual), "infinity bound"
    print(f"TORSION_NATIVE_OK: {len(actual)} complete point sets and invariant factors")


def test_unsupported_cells():
    result = invoke("2 0\n3 0\n2 1\n3 1\n")
    assert result.returncode == 0, result.stderr
    for row in csv.DictReader(io.StringIO(result.stdout)):
        assert row["status"] == "unsupported", row
        assert all(row[field] == "" for field in ("torsion", "torsion_points", "torsion_orders")), row


if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
    test_unsupported_cells()
