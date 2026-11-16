import csv
import ctypes
import ctypes.util
import io
from pathlib import Path
import subprocess
import sys
from _support import BINARY, ROOT, fake_magma, invoke


def decompress(data):
    # Streaming frames do not carry an uncompressed size.
    library = ctypes.CDLL(ctypes.util.find_library("zstd"))
    library.ZSTD_decompressBound.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    library.ZSTD_decompressBound.restype = ctypes.c_ulonglong
    library.ZSTD_decompress.argtypes = [ctypes.c_void_p, ctypes.c_size_t,
                                       ctypes.c_void_p, ctypes.c_size_t]
    library.ZSTD_decompress.restype = ctypes.c_size_t
    library.ZSTD_isError.argtypes = [ctypes.c_size_t]
    bound = library.ZSTD_decompressBound(data, len(data))
    assert bound < 100_000_000, "invalid or excessive zstd frame"
    target = ctypes.create_string_buffer(bound)
    size = library.ZSTD_decompress(target, bound, data, len(data))
    assert not library.ZSTD_isError(size), "zstd decompression failed"
    return target.raw[:size].decode("ascii")


def test_equations(reference):
    expected_text = Path(reference).read_text()
    expected = list(csv.DictReader(io.StringIO(expected_text)))
    text = "".join(f"{r['q']} {r['key']}\n" for r in expected)
    with fake_magma() as (env, marker):
        result = invoke(text, diagnostic=False, env=env)
        assert not marker.exists(), "native command started Magma"
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
    print(f"RECORDS_NATIVE_OK: {len(rows)} rows equal Magma byte for byte")


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
    print('RECORDS_ERRORS_OK: input and output failures return nonzero; preceding rows remain decodable')


if __name__ == '__main__':
    test_equations(sys.argv[1])
    test_execution_errors()
