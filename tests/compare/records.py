import csv
import io
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from _support import BINARY, ROOT, invoke


def decompress(data):
    with tempfile.NamedTemporaryFile(dir=ROOT / "build", suffix=".zst") as frame:
        frame.write(data)
        frame.flush()
        listing = subprocess.run(["zstd", "--list", frame.name], check=True,
                                 stdout=subprocess.PIPE, text=True).stdout
        fields = listing.splitlines()[1].split()
        assert fields[:2] == ["1", "0"], "expected one zstd frame without skips"
        assert "XXH64" in fields, "zstd frame lacks a checksum"
        decoded = subprocess.run(["zstd", "-dc", frame.name], check=True,
                                 stdout=subprocess.PIPE).stdout.decode("ascii")
    assert decoded.endswith("\n"), "CSV lacks a final newline"
    return decoded


def test_equations(reference):
    expected_text = Path(reference).read_text()
    expected = list(csv.DictReader(io.StringIO(expected_text)))
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    result = invoke(text, diagnostic=False)
    assert result.returncode == 0, result.stderr
    decoded = decompress(result.stdout)
    assert decoded == expected_text, "native CSV differs from Magma"
    rows = list(csv.DictReader(io.StringIO(decoded)))
    assert len({(r['kind'], r['minimality']) for r in rows}) == 6
    assert any(r['admissible_sha'] == '1 4 9' for r in rows)
    assert any(r['sha_reg'] == '3/4' for r in rows)
    assert {r['refinement_status'] for r in rows} >= {'', 'unsupported', 'not_needed'}
    assert any(r['torsion'] == '[]' for r in rows)
    assert any(r['torsion'] == '' for r in rows)


def test_execution_errors():
    invalid = invoke('2 0\n2 bad\n3 0\n', diagnostic=False)
    assert invalid.returncode != 0 and b'line 2' in invalid.stderr, invalid.stderr
    rows = list(csv.DictReader(io.StringIO(decompress(invalid.stdout))))
    assert len(rows) == 1 and rows[0]['q'] == '2'
    empty = invoke(' \n\n', diagnostic=False)
    assert empty.returncode == 0
    assert not list(csv.DictReader(io.StringIO(decompress(empty.stdout))))
    with open('/dev/full', 'wb') as sink:
        failure = subprocess.run([BINARY, 'invariants'], input=b'2 0\n',
                                 stdout=sink, stderr=subprocess.PIPE, cwd=ROOT)
    assert failure.returncode != 0 and b'write' in failure.stderr, failure.stderr


if __name__ == '__main__':
    test_equations(sys.argv[1])
    test_execution_errors()
