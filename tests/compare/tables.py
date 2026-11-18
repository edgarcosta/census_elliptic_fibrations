import csv
import hashlib
import io
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"


def count(rows):
    return subprocess.run([BINARY, "count"], input=rows, text=True, capture_output=True)


def test_all_entries(referencepath):
    lines = Path(referencepath).read_text().splitlines()
    assert len(lines) == 72, len(lines)
    hashes = dict(
        reversed(line.split()) for line in
        (ROOT / "tests/fixtures/trace_tables/tables.sha256").read_text().splitlines()
    )
    total = 0
    for offset in range(0, len(lines), 3):
        p, degree = map(int, lines[offset].split())
        expected = [1, *map(int, lines[offset + 2].split())]
        parameters = ["zero", *map(str, range(p**degree - 1))]
        request = "".join(f"{p} {degree} {e}\n" for e in parameters)
        result = count(request)
        assert result.returncode == 0, result.stderr
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        assert len(rows) == p**degree, (p, degree, len(rows))
        actual = []
        for parameter, row in zip(parameters, rows):
            assert (int(row["p"]), int(row["degree"]), row["parameter"]) == (
                p, degree, parameter
            ), row
            trace = int(row["trace"])
            assert int(row["points"]) == p**degree + 1 - trace, row
            actual.append(trace)
        assert actual == expected, (p, degree)
        filename = f"f{p}-d{degree}.txt"
        canonical = "\n".join(lines[offset:offset + 3]) + "\n"
        assert hashlib.sha256(canonical.encode()).hexdigest() == hashes[filename]
        total += len(expected) - 1
    assert total == 805326, total


def test_rejected_inputs():
    for row in ("5 1 0", "2 0 0", "3 13 0", "2 1 1", "3 1 2",
                "3 12 531440", "2 1 -1", "2 1 18446744073709551616",
                "2 1", "2 1 0 extra", "2 x 0"):
        result = count(row + "\n")
        assert result.returncode != 0, row
        assert result.stderr, row
        assert len(result.stdout.splitlines()) <= 1, (row, result.stdout)


def test_in_memory_tables():
    result = count("3 1 zero\n2 1 0\n3 1 1\n3 2 1\n3 1 zero\n")
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(io.StringIO(result.stdout)))
    assert [int(row["trace"]) for row in rows] == [1, -1, 1, 4, 1], rows
    assert [int(row["points"]) for row in rows] == [3, 4, 3, 6, 3], rows


if __name__ == "__main__":
    test_all_entries(sys.argv[1])
    test_rejected_inputs()
    test_in_memory_tables()
