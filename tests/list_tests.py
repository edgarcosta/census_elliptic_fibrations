#!/usr/bin/env python3
import collections
from pathlib import Path
import sys


paths = sorted(
    path.as_posix()
    for directory in ("tests", "tests/native", "examples")
    for path in Path(directory).glob("*.m")
    if path.is_file() and not path.name.startswith("_")
    and path != Path("tests/run_tests.m")
)
counts = collections.Counter(Path(path).name for path in paths)
duplicates = sorted(name for name, count in counts.items() if count > 1)
if duplicates:
    sys.exit("Duplicate test basenames: " + ", ".join(duplicates))
print("\n".join(paths))
