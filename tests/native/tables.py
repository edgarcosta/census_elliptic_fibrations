import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"


def count(rows, cache=None):
    arguments = [BINARY, "count"]
    if cache is not None:
        arguments += ["--table-cache", cache]
    return subprocess.run(arguments, input=rows, text=True, capture_output=True)


def direct_traces(p, polynomial):
    degree = len(polynomial) - 1
    size = p**degree

    def digits(a):
        result = []
        for _ in range(degree):
            result.append(a % p)
            a //= p
        return result

    def encode(values):
        return sum((v % p) * p**i for i, v in enumerate(values))

    def add(a, b):
        return encode([x + y for x, y in zip(digits(a), digits(b))])

    def multiply(a, b):
        product = [0] * (2 * degree - 1)
        for i, x in enumerate(digits(a)):
            for j, y in enumerate(digits(b)):
                product[i + j] += x * y
        for i in range(len(product) - 1, degree - 1, -1):
            for j in range(degree):
                product[i - degree + j] -= product[i] * polynomial[j]
        return encode(product[:degree])

    squares = [multiply(x, x) for x in range(size)]
    cubes = [multiply(x, squares[x]) for x in range(size)]
    # Enumerate both affine coordinates, independently of character transforms.
    traces = []
    for a in range(size):
        affine = sum(
            add(squares[y], multiply(x, y)) == add(cubes[x], a)
            if p == 2 else squares[y] == add(add(cubes[x], squares[x]), a)
            for x in range(size) for y in range(size)
        )
        traces.append(size - affine)
    alpha = p if degree > 1 else (-polynomial[0]) % p
    powers, a = [], 1
    for _ in range(size - 1):
        powers.append(traces[a])
        a = multiply(a, alpha)
    assert a == 1
    return [traces[0], *powers]


def test_all_entries(referencepath):
    lines = Path(referencepath).read_text().splitlines()
    assert len(lines) == 72, len(lines)
    hashes = dict(
        reversed(line.split()) for line in
        (ROOT / "tests/native/data/trace_tables/tables.sha256").read_text().splitlines()
    )
    total = 0
    with tempfile.TemporaryDirectory(dir=ROOT / "artifacts") as directory:
        cache = Path(directory)
        for offset in range(0, len(lines), 3):
            p, degree = map(int, lines[offset].split())
            polynomial = list(map(int, lines[offset + 1].split()))
            expected = [1, *map(int, lines[offset + 2].split())]
            parameters = ["zero", *map(str, range(p**degree - 1))]
            request = "".join(f"{p} {degree} {e}\n" for e in parameters)
            result = count(request, cache)
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
            if degree <= 3:
                assert actual == direct_traces(p, polynomial), (p, degree)
            filename = f"f{p}-d{degree}.txt"
            canonical = "\n".join(lines[offset:offset + 3]) + "\n"
            assert (cache / filename).read_text() == canonical, filename
            assert hashlib.sha256(canonical.encode()).hexdigest() == hashes[filename]
            # A fresh process certifies cached data before using it to count.
            cached = count(f"{p} {degree} zero\n{p} {degree} 0\n", cache)
            assert cached.returncode == 0, cached.stderr
            assert cached.stdout.splitlines() == result.stdout.splitlines()[:3]
            total += len(expected) - 1
        for p, degree in ((2, 3), (3, 12)):
            path = cache / f"f{p}-d{degree}.txt"
            original = path.read_text()
            parts = original.splitlines()
            entries = parts[2].split()
            entries[len(entries) // 2] = str(int(entries[len(entries) // 2]) + 1)
            parts[2] = " ".join(entries)
            path.write_text("\n".join(parts) + "\n")
            result = count(f"{p} {degree} 0\n", cache)
            assert result.returncode != 0, (p, degree)
            assert "certificate" in result.stderr, result.stderr
            assert len(result.stdout.splitlines()) <= 1, result.stdout
            assert path.read_text() == "\n".join(parts) + "\n"
            path.write_text(original)
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
    result = count("3 1 zero\n2 1 0\n3 1 1\n3 1 zero\n")
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(io.StringIO(result.stdout)))
    assert [int(row["trace"]) for row in rows] == [1, -1, 1, 1], rows
    assert [int(row["points"]) for row in rows] == [3, 4, 3, 3], rows


def test_concurrent_cold_cache():
    with tempfile.TemporaryDirectory(dir=ROOT / "artifacts") as directory:
        with ThreadPoolExecutor(max_workers=4) as workers:
            for round_number in range(20):
                cache = Path(directory) / str(round_number)
                results = list(workers.map(lambda _: count("3 12 0\n", cache), range(4)))
                for result in results:
                    assert result.returncode == 0, (round_number, result.stderr)
                    assert result.stdout == (
                        "p,degree,parameter,trace,points\n3,12,0,658,530784\n"
                    ), (round_number, result.stdout)
                warm = count("3 12 0\n", cache)
                assert warm.returncode == 0, (round_number, warm.stderr)
                assert warm.stdout == results[0].stdout, (round_number, warm.stdout)


if __name__ == "__main__":
    test_all_entries(sys.argv[1])
    test_rejected_inputs()
    test_in_memory_tables()
    test_concurrent_cold_cache()
