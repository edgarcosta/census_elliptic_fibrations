import csv
import io
import sys
from _support import fake_magma, invoke


LOCAL_FIELDS = {"fibres", "twist_fibres", "degree", "sign", "height",
                "twist_height", "tamagawa_exponent"}


def fibres(text):
    if text == "[]":
        return []
    rows = []
    for entry in text.split(";"):
        place, degree, kind, index, vd, conductor, c, trace, cycles = entry.split(":")
        rows.append((place, int(degree), kind, int(index), int(vd),
                     int(conductor), int(c), int(trace),
                     [] if cycles == "[]" else list(map(int, cycles.split()))))
    return rows


def check_identities(text, height):
    rows = fibres(text)
    assert sum(r[1] * r[4] for r in rows) == 12 * height, rows
    for place, d, kind, n, vd, conductor, c, trace, cycles in rows:
        components = n if kind == "I" else n + 5 if kind == "I*" else {
            "II": 1, "III": 2, "IV": 3, "IV*": 7, "III*": 8, "II*": 9,
        }[kind]
        assert conductor == vd + 1 - components, (place, kind, conductor, vd)
        assert sum(cycles) == d * (components - 1), (place, kind, cycles)
        assert cycles == sorted(cycles) and all(k % d == 0 for k in cycles)


def test_magma_comparison(reference):
    with open(reference, newline="") as source:
        expected = list(csv.DictReader(source))
    assert expected
    for p in ("2", "3"):
        sample = [r for r in expected if r["q"] == p]
        keys = {r["key"] for r in sample}
        assert len(keys) >= 100, (p, "distinct keys", len(keys))
        reductions = [f for r in sample for field in ("fibres", "twist_fibres")
                      for f in fibres(r[field])]
        families = {f[2] for f in reductions}
        assert families == {
            "I", "I*", "II", "III", "IV", "IV*", "III*", "II*",
        }, (p, families)
        assert any(f[2] == "I*" and f[3] > 0 for f in reductions), p
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    # A runtime Magma fallback must fail even when the test host has Magma.
    with fake_magma() as (env, marker):
        result = invoke(text, env=env)
        assert not marker.exists(), "native command started Magma"
    assert result.returncode == 0, result.stderr
    actual = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(actual) == len(expected), (len(actual), len(expected))
    for got, want in zip(actual, expected):
        context = f"q={want['q']}, key={want['key']}"
        for field, value in want.items():
            assert got[field] == value, (context, field, got.get(field), value)
        assert got["status"] == "partial", context
        assert int(got["height"]) == 2, context
        check_identities(got["fibres"], 2)
        check_identities(got["twist_fibres"], int(got["twist_height"]))
    assert {r["sign"] for r in actual} == {"1", "-1"}
    assert {r["twist_height"] for r in actual} >= {"1", "2", "4"}


def test_unsupported_cells():
    result = invoke("2 0\n3 0\n2 1\n3 1\n")
    assert result.returncode == 0, result.stderr
    for row in csv.DictReader(io.StringIO(result.stdout)):
        assert row["status"] == "unsupported", row
        assert all(row[field] == "" for field in LOCAL_FIELDS), row


if __name__ == "__main__":
    test_magma_comparison(sys.argv[1])
    test_unsupported_cells()
