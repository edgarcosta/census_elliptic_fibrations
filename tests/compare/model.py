import csv
import io
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"


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
        supported = want["kind"] == "smooth" and want["minimality"] == "everywhere minimal"
        assert got["status"] == ("partial" if supported else "unsupported"), context
        assert set(got) == set(want) | {"status"}, got
    for p in ("2", "3"):
        cells = {(row["kind"], row["minimality"]) for row in actual if row["q"] == p}
        assert len(cells) == 6, cells


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
    result = invoke("2 0\n", "--unknown")
    assert result.returncode != 0 and result.stderr, result



if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
    test_invalid_input()
