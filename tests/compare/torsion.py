import csv
import io
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from _support import invoke


def test_magma_comparison(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert len(expected) == 34, len(expected)
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    result = invoke(text)
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected)
    for got, want in zip(actual, expected):
        context = (want["q"], want["key"])
        assert got["torsion"] == want["torsion"], context
        points = got["torsion_points"].split(";")
        wanted = want["torsion_points"].split(";")
        assert len(points) == len(set(points)) == len(wanted), context
        assert set(points) == set(wanted), context
        assert points[0] == "O", context
    assert {r["torsion"] for r in actual} >= {"[]", "2", "3", "6", "2 4"}
    assert {r["torsion"] for r in actual if r["q"] == "3"} >= {"4", "2 2", "2 6"}
    assert any("1:4:1:8" in r["torsion_points"].split(";") for r in actual), "degree-one denominator"
    assert any("64:1:512:1" in r["torsion_points"].split(";") for r in actual), "infinity bound"


if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
