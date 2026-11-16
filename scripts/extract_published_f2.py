"""Extract 16 pinned rows from a local F_2 census table checkout.

From the repository root:
    python3 scripts/extract_published_f2.py TABLE_CHECKOUT OUTPUT.csv
TABLE_CHECKOUT must contain the pinned revision; OUTPUT.csv is the sample path.
See tests/data/published-f2/README.md for fixture provenance.
"""

import csv
import io
from pathlib import Path
import subprocess
import sys
import zipfile


REVISION = "9d1b3b85b702811da14d0e53a3d768356e998f27"
ARCHIVE_BLOB = "0f1839c0e0f8746ae8a93b9fda456c250aebcb16"
CLASS_SIZES = {"0": 173977, "s^2": 524493, "s*t": 510215,
               "s^2 + s*t + t^2": 40884}
FIELDS = ["key", "a1", "mw_torsion", "singular_fibres", "L_function",
          "L_sign", "analytic_rank"]


def extract(checkout, output):
    blob = subprocess.check_output(
        ["git", "-C", checkout, "rev-parse", f"{REVISION}:census.csv.zip"],
        text=True).strip()
    if blob != ARCHIVE_BLOB:
        raise ValueError("census.csv.zip does not match the pinned Git blob")
    archive = subprocess.check_output(
        ["git", "-C", checkout, "show", f"{REVISION}:census.csv.zip"])
    positions = {a: {0, n // 3, (2 * n) // 3, n - 1}
                 for a, n in CLASS_SIZES.items()}
    counts = dict.fromkeys(CLASS_SIZES, 0)
    sample = []
    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        with zipped.open("census.csv") as source:
            for index, row in enumerate(csv.DictReader(io.TextIOWrapper(source)), 1):
                a = row["a1"]
                if counts[a] in positions[a]:
                    sample.append({"source_row": index,
                                   **{field: row[field] for field in FIELDS}})
                counts[a] += 1
    if counts != CLASS_SIZES or len(sample) != 16:
        raise ValueError("The census.csv class sizes or sample size differ from the pinned values")
    with Path(output).open("w", newline="", encoding="ascii") as target:
        writer = csv.DictWriter(target, fieldnames=["source_row", *FIELDS],
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(sample)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: extract_published_f2.py TABLE_CHECKOUT OUTPUT.csv")
    extract(sys.argv[1], sys.argv[2])
