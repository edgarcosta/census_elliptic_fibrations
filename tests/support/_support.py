from contextlib import contextmanager
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("CENSUS_NATIVE", ROOT / "artifacts/native/census-native"))


@contextmanager
def fake_magma():
    with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2] / "artifacts") as directory:
        sentinel = Path(directory) / "magma"
        marker = Path(directory) / "called"
        sentinel.write_text(f"#!/bin/sh\n: > '{marker}'\nexit 99\n")
        sentinel.chmod(0o755)
        yield {**os.environ, "PATH": directory, "MAGMA": str(sentinel)}, marker


def invoke(text, *options, **kwargs):
    return subprocess.run([BINARY, "invariants", "--diagnostic", *options],
                          input=text, capture_output=True, text=True, cwd=ROOT, **kwargs)
