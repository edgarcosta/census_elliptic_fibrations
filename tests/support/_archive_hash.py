"""Hash canonical orbit records for comparison with the test constants.

Records are (I, smooth, minimal), with binary flags and minimal=0 for singular
models. The disjoint buckets are smooth_minimal, smooth_nonminimal and
delta_zero; all is their union. Sort full tuples, preserving multiplicity.

Each prefix digest hashes this ASCII stream with LF endings:
    census-orbit-v1|q|canonical_prefix_integer|bucket|record_count
    coefficient_digits_in_fixed_width|smooth_0_or_1|minimal_0_or_1
    ...
Coefficient width is 32 for F_2 and 27 for F_3, including leading zeroes.
An empty bucket hashes its header with count zero. Combine prefix digests
in increasing numerical prefix order:
    census-multiset-v1|q|bucket
    canonical_prefix_integer|record_count|prefix_bucket_digest
    ...
Counts and uniqueness are checked separately. Serialization and input order
do not affect the comparison; SHA-256 collisions remain its hash limit.
"""

import hashlib
from pathlib import Path
import struct
import sys


BUCKETS = ("all", "smooth_minimal", "smooth_nonminimal", "delta_zero")


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


def check_prefix(q, prefix, path):
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
        print(f"{prefix}\t{bucket}\t{actual[0]}\t{actual[1]}")


def check_combined(q):
    actual = {}
    for line in sys.stdin:
        prefix, bucket, count, sha = line.rstrip("\n").split("\t")
        key = q, int(prefix), bucket
        assert key not in actual
        actual[key] = int(count), sha
    labels = {key[1] for key in actual}
    assert len(labels) == (88 if q == 2 else 2)
    assert set(actual) == {(q, label, bucket) for label in labels for bucket in BUCKETS}
    for bucket in BUCKETS:
        entries = [(key[1], *value) for key, value in actual.items() if key[2] == bucket]
        count, sha = combined(q, bucket, entries)
        print(f"{bucket}\t{count}\t{sha}")


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
