from contextlib import contextmanager
import os
from pathlib import Path
import tempfile


@contextmanager
def fake_magma():
    with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2] / "artifacts") as directory:
        sentinel = Path(directory) / "magma"
        marker = Path(directory) / "called"
        sentinel.write_text(f"#!/bin/sh\n: > '{marker}'\nexit 99\n")
        sentinel.chmod(0o755)
        yield {**os.environ, "PATH": directory, "MAGMA": str(sentinel)}, marker
