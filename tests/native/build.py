from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]


def test_build_and_usage():
    build = subprocess.run(
        ["make", "census-native"], cwd=ROOT, capture_output=True, text=True
    )
    assert build.returncode == 0, build.stdout + build.stderr
    binary = ROOT / "artifacts/native/census-native"
    for option in ("--help", "-h", "help"):
        result = subprocess.run([binary, option], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert "Usage: census-native" in result.stdout, result.stdout
        assert "help" in result.stdout, result.stdout
    for arguments in ([], ["unknown"], ["--help", "extra"]):
        result = subprocess.run([binary, *arguments], capture_output=True, text=True)
        assert result.returncode != 0, arguments
        assert "Usage: census-native" in result.stderr, result.stderr


if __name__ == "__main__":
    test_build_and_usage()
