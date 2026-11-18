from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"


def invoke(text, *options):
    return subprocess.run([BINARY, "invariants", "--diagnostic", *options],
                          input=text, capture_output=True, text=True, cwd=ROOT)
