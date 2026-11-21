from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "artifacts/native/census-native"


def invoke(text, *options, diagnostic=True):
    mode = ["--diagnostic"] if diagnostic else []
    return subprocess.run([BINARY, "invariants", *mode, *options],
                          input=text if diagnostic else text.encode("ascii"),
                          capture_output=True, text=diagnostic, cwd=ROOT)
