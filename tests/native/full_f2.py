import collections
import csv
from fractions import Fraction
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import sqlite3
import struct
import subprocess
import time
import zipfile

from _support import ROOT, fake_magma, invoke
from candidates import check_candidates
from derived import check_identities as check_derived, integers
from local import check_identities as check_local
from records import decompress


DATA = ROOT / "tests/native/data/full_f2"
OUT = ROOT / "artifacts/full-f2"
A1_CLASS = {0: "zero", 1: "square", 4: "square", 5: "square",
            2: "split", 3: "split", 6: "split", 7: "irreducible"}
SHARED = ("a1", "mw_torsion", "singular_fibres", "L_function", "L_sign",
          "analytic_rank", "P2", "zeta_sign", "rho_F2", "rho_geom",
          "artin_height", "br_disc", "pts_F2", "pts_F4", "pts_F8",
          "pts_F16", "pts_F32")


def digest(path):
    with open(path, "rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def published_rows():
    with zipfile.ZipFile(DATA / "census.csv.zip") as archive:
        with archive.open("census.csv") as source:
            yield from csv.DictReader(io.TextIOWrapper(source, encoding="ascii"))


def published_values(row):
    mask = 0
    if row["a1"] != "0":
        for term in row["a1"].replace(" ", "").split("+"):
            mask ^= {"s^2": 4, "s*t": 2, "t^2": 1}[term]
    assert mask == int(row["key"]) >> 29, row
    values = dict(row)
    values["a1"] = A1_CLASS[mask]
    values["mw_torsion"] = ([] if row["mw_torsion"] == "trivial" else
                            [int(v.strip()[2:]) for v in row["mw_torsion"].split("x")])
    values["singular_fibres"] = sorted(
        [int(d), int(f), int(c), kind]
        for d, f, c, kind in (s.split(":") for s in row["singular_fibres"].split(";") if s))
    for field in ("L_function", "P2"):
        values[field] = integers(row[field])
    for field in SHARED[4:]:
        if field != "P2" and values[field] != "inf":
            values[field] = int(values[field])
    return [values[field] for field in SHARED]


def native_values(row):
    fibres = [f.split(":") for f in row["fibres"].split(";")]
    values = [A1_CLASS[int(row["key"]) >> 29], integers(row["torsion"]),
              sorted([int(f[1]), int(f[5]), int(f[6]),
                      "I" + f[3] + ("*" if f[2] == "I*" else "")
                      if f[2] in ("I", "I*") else f[2]]
                     for f in fibres),
              integers(row["L_polynomial"]), int(row["sign"]), int(row["rank"]),
              integers(row["P2"]), int(row["zeta_sign"]),
              int(row["picard_number"]), int(row["geometric_picard_number"]),
              "inf" if row["brauer_height"] == "infinity" else int(row["brauer_height"]),
              int(row["artin_tate_product"])]
    return values + integers(row["point_counts"])


def check_row(record, diagnostic, key):
    context = f"F_2 key {key}"
    assert (record["q"], record["key"], record["schema_version"], record["status"]) == (
        "2", str(key), "1", "ok"), (context, record)
    assert record["kind"] == diagnostic["kind"] == "smooth", context
    assert record["minimality"] == diagnostic["minimality"] == "everywhere minimal", context
    for field in ("degree", "sign", "lower_L", "torsion",
                  "admissible_sha", "split_method"):
        assert record[field] == diagnostic[field], (context, field)
    assert record["fibres"] == ";".join(
        ":".join(f.split(":")[:8]) for f in diagnostic["fibres"].split(";")), context
    assert record["message"] == record["d2"] == record["sha2_dimensions"] == "", context
    assert record["refinement_status"] == (
        "not_needed" if record["split_method"] == "constraints" else "unsupported"), context
    coefficients = integers(diagnostic["L_polynomial"])
    degree, sign = int(record["degree"]), int(record["sign"])
    assert len(coefficients) == degree + 1 and coefficients[0] == 1, context
    assert sign in (-1, 1), context
    assert integers(record["lower_L"]) == coefficients[1:1 + degree//2], context
    assert all(coefficients[degree-i] == sign*2**(degree-2*i)*coefficients[i]
               for i in range(degree//2 + 1)), context
    assert (sign == -1) == (int(diagnostic["rank"]) % 2 == 1), context
    assert degree == sum(int(f.split(":")[1])*int(f.split(":")[5])
                         for f in record["fibres"].split(";")) - 4, context
    codes = [sum(c*2**i for i, c in enumerate(integers(diagnostic[a])))
             for a in ("a1", "a2", "a3", "a4", "a6")]
    assert codes == [key >> 29, 0, key >> 22 & 127, key >> 13 & 511, key & 8191], context
    fibres = [f.split(":") for f in diagnostic["fibres"].split(";")]
    cycles = [n for f in fibres for n in integers(f[8])]
    p2 = coefficients[:]
    for n in [1, 1] + cycles:
        product = p2 + [0]*n
        for i, c in enumerate(p2):
            product[i+n] -= 2**n*c
        p2 = product
    assert p2 == integers(diagnostic["P2"]), context
    quotient = coefficients[:]
    for _ in range(int(diagnostic["rank"])):
        divided = [quotient[0]]
        for c in quotient[1:]:
            divided.append(c + 2*divided[-1])
        assert divided.pop() == 0, context
        quotient = divided
    d = len(quotient)-1
    leading = Fraction(2*sum(c*2**(d-i) for i, c in enumerate(quotient)), 2**d)
    torsion = math.prod(integers(record["torsion"]))
    assert Fraction(record["sha_reg"]) == leading*torsion**2 / math.prod(int(f[6]) for f in fibres), context
    check_local(diagnostic["fibres"], 2)
    check_derived(diagnostic)
    check_candidates(diagnostic, diagnostic)
    assert integers(record["admissible_sha"]), context
    assert Fraction(record["sha_reg"]) == Fraction(diagnostic["sha_reg"]) > 0, context


def check_f3(reference):
    expected = reference.read_text(encoding="ascii")
    rows = list(csv.DictReader(io.StringIO(expected)))
    assert len(rows) == 5 and all(r["q"] == "3" and r["status"] == "ok" for r in rows)
    with fake_magma() as (env, marker):
        result = invoke("".join(f"3 {r['key']}\n" for r in rows), diagnostic=False, env=env)
        assert not marker.exists(), "native command started Magma"
    assert result.returncode == 0, result.stderr
    assert decompress(result.stdout) == expected, "F_3 records differ from Magma"
    return len(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    manifest = json.loads((DATA / "sources.json").read_text())
    assert digest(DATA / "census.csv.zip") == manifest["archive_sha256"], "census archive hash; run git lfs pull"
    with zipfile.ZipFile(DATA / "census.csv.zip") as archive:
        with archive.open("census.csv") as source:
            assert hashlib.file_digest(source, "sha256").hexdigest() == manifest["csv_sha256"]
    keypath, reference = OUT / "representatives.csv", OUT / "f3-reference.csv"
    with (OUT / "magma.log").open("w") as log:
        subprocess.run([os.environ.get("MAGMA", "magma"), "-bn", f"keypath:={keypath}",
                        f"referencepath:={reference}", "tests/native/full_f2_reference.m"],
                       stdout=log, stderr=subprocess.STDOUT, check=True, cwd=ROOT)
    keys = []
    with keypath.open() as source:
        for row in csv.DictReader(source):
            if row["cell"] == "0":
                keys.append(int(row["key"]))
    keys.sort()
    assert len(keys) == len(set(keys)) == manifest["native_rows"]
    population_hash = hashlib.sha256(b"".join(struct.pack("<I", k) for k in keys)).hexdigest()
    assert population_hash == manifest["population_sha256"], "F_2 orbit population changed"
    f3_count = check_f3(reference)
    with (OUT / "published-keys.txt").open("w") as target:
        count = 0
        for row in published_rows():
            target.write(row["key"] + "\n")
            count += 1
    assert count == manifest["published_rows"]
    with (OUT / "orbit.log").open("w") as log:
        subprocess.run([ROOT / "artifacts/native/f2-orbit-match", keypath,
                        OUT / "published-keys.txt", OUT / "matches.csv", OUT / "hits.bin"],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    correspondence = {}
    with (OUT / "matches.csv").open() as source, (OUT / "unmatched-published.csv").open("w") as output:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(output, fieldnames=reader.fieldnames)
        writer.writeheader()
        unmatched_published = 0
        for row in reader:
            if row["n_smooth_min"] != "1" or row["other_cells"]:
                writer.writerow(row)
                unmatched_published += 1
            else:
                published, native = int(row["published_key"]), int(row["native_key"])
                assert published not in correspondence, ("duplicate published key", published)
                correspondence[published] = native
    assert len(set(correspondence.values())) == len(correspondence), "noninjective orbit match"
    database = OUT / "comparison.sqlite"
    database.unlink(missing_ok=True)
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE expected (key INTEGER PRIMARY KEY, published INTEGER, fields TEXT)")
        connection.executemany("INSERT INTO expected VALUES (?, ?, ?)",
                               ((correspondence[int(row["key"])], int(row["key"]),
                                 json.dumps(published_values(row)))
                                for row in published_rows() if int(row["key"]) in correspondence))
        connection.commit()
        del correspondence
        summary = {"archive_sha256": manifest["archive_sha256"], "population_sha256": population_hash,
                   "published_rows": count, "native_rows": 0, "identity_rows": 0,
                   "matched_rows": 0, "unmatched_published": unmatched_published,
                   "unmatched_native": 0, "f3_controls": f3_count,
                   "agree_by_column": collections.Counter(), "disagree_by_column": collections.Counter(),
                   "machine": platform.node(), "workers": 1}
        with (OUT / "unmatched-native.csv").open("w") as missing, (OUT / "disagreements.csv").open("w") as bad:
            missing_writer, bad_writer = csv.writer(missing), csv.writer(bad)
            missing_writer.writerow(["key"])
            bad_writer.writerow(["published_key", "native_key", "field", "published", "native"])
            with fake_magma() as (env, marker):
                for begin in range(0, len(keys), 10000):
                    batch = keys[begin:begin+10000]
                    text = "".join(f"2 {key}\n" for key in batch)
                    records = invoke(text, diagnostic=False, env=env)
                    diagnostic = invoke(text, env=env)
                    assert records.returncode == 0, records.stderr
                    assert diagnostic.returncode == 0, diagnostic.stderr
                    compact_rows = list(csv.DictReader(io.StringIO(decompress(records.stdout))))
                    diagnostic_rows = list(csv.DictReader(io.StringIO(diagnostic.stdout)))
                    assert len(compact_rows) == len(diagnostic_rows) == len(batch), begin
                    for key, record, expanded in zip(batch, compact_rows, diagnostic_rows):
                        assert expanded["key"] == str(key) and expanded["q"] == "2", key
                        check_row(record, expanded, key)
                        summary["native_rows"] += 1
                        summary["identity_rows"] += 1
                        expected = connection.execute("SELECT published, fields FROM expected WHERE key=?", (key,)).fetchone()
                        if expected is None:
                            missing_writer.writerow([key])
                            summary["unmatched_native"] += 1
                            continue
                        summary["matched_rows"] += 1
                        actual = native_values(expanded)
                        assert len(actual) == len(SHARED)
                        for field, want, got in zip(SHARED, json.loads(expected[1]), actual):
                            if want == got:
                                summary["agree_by_column"][field] += 1
                            else:
                                summary["disagree_by_column"][field] += 1
                                bad_writer.writerow([expected[0], key, field, want, got])
                    print(f"checked {summary['native_rows']}/{len(keys)} F_2 rows", flush=True)
                assert not marker.exists(), "native command started Magma"
        summary["unique_completed_rows"] = summary["native_rows"]
        summary["fields_compared"] = sum(summary["agree_by_column"].values()) + sum(summary["disagree_by_column"].values())
        summary["wall_seconds"] = time.monotonic() - start
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        assert summary["native_rows"] == summary["identity_rows"] == manifest["native_rows"]
        assert summary["matched_rows"] == count and unmatched_published == 0
        assert summary["unmatched_native"] == manifest["expected_unmatched_native"]
        assert not summary["disagree_by_column"], summary
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
