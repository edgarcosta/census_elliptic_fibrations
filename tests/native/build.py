from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def test_native_test_dependencies():
    (ROOT / "artifacts").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=ROOT / "artifacts") as directory:
        root = Path(directory)
        shutil.copyfile(ROOT / "Makefile", root / "Makefile")
        (root / "native").mkdir()
        (root / "native/main.cpp").touch()
        (root / "tests").mkdir()
        (root / "tests/native_probe.m").touch()
        for target in ("check/native_probe.m", "debug/native_probe.m"):
            result = subprocess.run(
                ["make", "--dry-run", target],
                cwd=root, capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            commands = result.stdout
            build = "-o artifacts/native/census-native"
            assert commands.count(build) == 1, commands
            invocation = "filename:=tests/native_probe.m"
            assert commands.index(build) < commands.index(invocation), commands


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
    test_native_test_dependencies()
    test_build_and_usage()
