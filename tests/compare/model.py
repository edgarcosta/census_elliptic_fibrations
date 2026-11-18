import csv
import io
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"
CHART = "a_i(t,1)"


def invoke(text, *options):
    return subprocess.run(
        [BINARY, "invariants", "--diagnostic", *options], input=text,
        capture_output=True, text=True, cwd=ROOT,
    )


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
        assert got["chart"] == CHART, context
        supported = want["kind"] == "smooth" and want["minimality"] == "everywhere minimal"
        assert got["status"] == ("partial" if supported else "unsupported"), context
        assert set(got) == set(want) | {"chart", "status"}, got
    for p in ("2", "3"):
        cells = {(row["kind"], row["minimality"]) for row in actual if row["q"] == p}
        assert len(cells) == 6, cells


def test_invalid_input():
    for p, bound in ((2, 2**32), (3, 3**27)):
        for key in ("-1", str(bound), str(bound + 1), str(2**64 - 1),
                    str(2**64), str(2**64 + 1), "9" * 200, "1x", "+1", "1.0"):
            result = invoke(f"{p} {key}\n")
            assert result.returncode != 0 and "key" in result.stderr.lower(), (p, key, result)
            assert not list(csv.DictReader(io.StringIO(result.stdout))), result.stdout
    for text in ("4 0\n", "9 0\n", "5 0\n", "-2 0\n", "2x 0\n",
                 "18446744073709551618 0\n", "2\n", "2 0 extra\n"):
        result = invoke(text)
        assert result.returncode != 0 and result.stderr, (text, result)
        assert not list(csv.DictReader(io.StringIO(result.stdout))), result.stdout
    for options in (("--chart", "a_i(1,t)"), ("--chart",), ("--unknown",)):
        result = invoke("2 0\n", *options)
        assert result.returncode != 0 and result.stderr, (options, result)
    result = invoke("3 0\n", "--chart", CHART)
    assert result.returncode == 0, result.stderr


if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
    test_invalid_input()
