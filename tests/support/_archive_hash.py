"""Compare canonical mathematical records with the archived multiset digests."""

import hashlib
from pathlib import Path
import struct
import sys


BUCKETS = ("all", "smooth_minimal", "smooth_nonminimal", "delta_zero")
FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/oldcode-archive/summaries.tsv"


def coefficients(value, q):
    digits = []
    for _ in range(32 if q == 2 else 27):
        digits.append(str(value % q))
        value //= q
    assert value == 0
    return "".join(reversed(digits))


def digest(q, prefix, bucket, records):
    header = f"census-orbit-v1|{q}|{prefix}|{bucket}|{len(records)}\n"
    result = hashlib.sha256(header.encode("ascii"))
    for value, smooth, minimal in sorted(records):
        row = f"{coefficients(value, q)}|{smooth}|{minimal}\n"
        result.update(row.encode("ascii"))
    return result.hexdigest()


def combined(q, bucket, entries):
    message = f"census-multiset-v1|{q}|{bucket}\n"
    message += "".join(f"{prefix}|{count}|{sha}\n"
                       for prefix, count, sha in sorted(entries))
    return sum(row[1] for row in entries), hashlib.sha256(message.encode("ascii")).hexdigest()


def fixture():
    lines = FIXTURE.read_text(encoding="ascii").splitlines()
    assert lines[0] == "q\told_prefix_index\tcanonical_prefix\tbucket\tcount\tsha256"
    assert lines[-1] == "END_SUMMARIES"
    boundary = lines.index("CELLS")
    prefixes, cells, ordinals = {}, {}, {}
    for line in lines[1:boundary]:
        q, ordinal, prefix, bucket, count, sha = line.split("\t")
        q, ordinal, prefix, count = map(int, (q, ordinal, prefix, count))
        assert q in (2, 3) and bucket in BUCKETS and count >= 0
        assert len(sha) == 64 and set(sha) <= set("0123456789abcdef")
        key = q, prefix, bucket
        assert key not in prefixes
        prefixes[key] = count, sha
        ordinals.setdefault((q, prefix), set()).add(ordinal)
    for line in lines[boundary + 1:-1]:
        q, bucket, count, sha = line.split("\t")
        key = int(q), bucket
        assert key not in cells
        cells[key] = int(count), sha
    for q, size in ((2, 88), (3, 2)):
        labels = {prefix for field, prefix in ordinals if field == q}
        assert len(labels) == size
        indices = [next(iter(ordinals[q, prefix])) for prefix in labels]
        assert all(len(ordinals[q, prefix]) == 1 for prefix in labels)
        assert set(indices) == (set(range(1, 89)) if q == 2 else {1, 934})
        assert q == 2 or labels == {0, 21870}
        for prefix in labels:
            assert prefixes[q, prefix, "all"][0] == sum(
                prefixes[q, prefix, bucket][0] for bucket in BUCKETS[1:])
        for bucket in BUCKETS:
            entries = [(prefix, *prefixes[q, prefix, bucket]) for prefix in labels]
            assert combined(q, bucket, entries) == cells[q, bucket]
    assert len(prefixes) == 360 and len(cells) == 8
    return prefixes, cells


def check_prefix(q, prefix, path):
    expected, _ = fixture()
    data = Path(path).read_bytes()
    records = list(struct.iter_unpack("<QBB", data))
    width = 22 if q == 2 else 13
    assert records and all(value // q**width == prefix for value, _, _ in records)
    assert len({value for value, _, _ in records}) == len(records), "duplicate orbit"
    groups = {bucket: [] for bucket in BUCKETS}
    groups["all"] = records
    for record in records:
        value, smooth, minimal = record
        assert 0 <= value < q**(32 if q == 2 else 27)
        assert smooth in (0, 1) and minimal in (0, 1) and (smooth or not minimal)
        bucket = "delta_zero" if not smooth else "smooth_minimal" if minimal else "smooth_nonminimal"
        groups[bucket].append(record)
    for bucket, group in groups.items():
        actual = len(group), digest(q, prefix, bucket, group)
        assert actual == expected[q, prefix, bucket], (q, prefix, bucket, actual, expected[q, prefix, bucket])
        print(f"{prefix}\t{bucket}\t{actual[0]}\t{actual[1]}")


def check_combined(q):
    expected, cells = fixture()
    actual = {}
    for line in sys.stdin:
        prefix, bucket, count, sha = line.rstrip("\n").split("\t")
        key = q, int(prefix), bucket
        assert key not in actual
        actual[key] = int(count), sha
    assert actual == {key: value for key, value in expected.items() if key[0] == q}
    for bucket in BUCKETS:
        entries = [(key[1], *value) for key, value in actual.items() if key[2] == bucket]
        assert combined(q, bucket, entries) == cells[q, bucket]
    print("ARCHIVE_OK")


if __name__ == "__main__":
    action, field, *args = sys.argv[1:]
    q = int(field)
    assert q in (2, 3)
    if action == "prefix":
        prefix, path = args
        check_prefix(q, int(prefix), path)
    else:
        assert action == "combined" and not args
        check_combined(q)
