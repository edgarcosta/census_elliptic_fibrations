import csv
import io
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from _support import invoke
from local import fibres
from derived import FIELDS as DERIVED_FIELDS
from candidates import FIELDS as CANDIDATE_FIELDS


def test_magma_comparison(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert expected
    result = invoke("".join(f"{row['q']} {row['key']}\n" for row in expected))
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected), (len(actual), len(expected))
    for got, want in zip(actual, expected):
        context = f"q={want['q']}, key={want['key']}"
        for field, value in want.items():
            assert got[field] == value, (context, field, got[field], value)
        supported = want["kind"] == "smooth" and want["minimality"] == "everywhere minimal"
        assert got["status"] == ("partial" if supported else "unsupported"), context
        assert set(got) == set(want) | set(DERIVED_FIELDS) | set(CANDIDATE_FIELDS) | {"torsion", "torsion_points", "lower_L", "L_polynomial", "trace_sums", "status"}, got
        if not supported:
            assert all(got[field] == "" for field in ("torsion", "torsion_points",
                       "lower_L", "L_polynomial", "trace_sums", *DERIVED_FIELDS, *CANDIDATE_FIELDS)), context
    for p in ("2", "3"):
        cells = {(row["kind"], row["minimality"]) for row in actual if row["q"] == p}
        assert len(cells) == 6, cells
    supported_rows = [row for row in actual if row["status"] == "partial"]
    assert {row["sign"] for row in supported_rows} == {"1", "-1"}
    assert {row["twist_height"] for row in supported_rows} >= {"1", "2", "4"}
    for p in ("2", "3"):
        reductions = [f for row in supported_rows if row["q"] == p
                      for field in ("fibres", "twist_fibres") for f in fibres(row[field])]
        assert {f[2] for f in reductions} == {"I", "I*", "II", "III", "IV", "IV*", "III*", "II*"}
        assert any(f[2] == "I*" and f[3] > 0 for f in reductions)


def test_invalid_input():
    for p, bound in ((2, 2**32), (3, 3**27)):
        for key in ("-1", str(bound), str(2**64), "1x"):
            result = invoke(f"{p} {key}\n")
            assert result.returncode != 0 and "key" in result.stderr.lower(), (p, key, result)
            assert not list(csv.DictReader(io.StringIO(result.stdout))), result.stdout
    for text in ("4 0\n", "2x 0\n", "4294967298 0\n", "2\n", "2 0 extra\n"):
        result = invoke(text)
        assert result.returncode != 0 and result.stderr, (text, result)
        assert not list(csv.DictReader(io.StringIO(result.stdout))), result.stdout
    for options in (("--unknown",), ("--backend",)):
        result = invoke("2 0\n", *options)
        assert result.returncode != 0 and result.stderr, (options, result)



if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
    test_invalid_input()
