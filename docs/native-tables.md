# Exact ordinary trace tables

`FibrationTraceTable(p,d)` constructs an exact table for p in [2,3] and extension
degree 1 <= d <= 12. Attach the package's `spec` before calling it.

The returned record retains:

- `Polynomial`: the defining monic polynomial over GF(p).
- `Field`: the extension defined by that polynomial.
- `Generator`: its primitive polynomial residue class alpha.
- `Ordinary`: the integer traces at alpha^e, stored in position e+1.
- `NodalTrace`: the trace 1 at parameter zero.

The equations are y^2+x*y=x^3+a in characteristic 2 and y^2=x^3+x^2+a in
characteristic 3. Nonzero parameters give smooth elliptic curves. At zero,
the point count is that of the singular projective cubic, with one point at
infinity. The table contains p^d-1 entries; the native prototype's two unused
trailing slots are omitted. In degree 1 over F_3, alpha is 2.

The F_3 implementation uses the [radix-3 certificate](native-mathematics.md):
M=1062907, omega=662797, inverse root omega^2=400109 and one final Q^-1
normalization, where Q=3^d. It centers every recovered coefficient and checks
its absolute value is at most Q. Since M>2Q, recovery is unique.
Construction takes O(d*p^d) arithmetic operations and O(p^d) storage.

Over F_2, let psi(z)=(-1)^Tr(z). For x != 0, division by x^2 turns the equation
into (y/x)^2+y/x=x+a/x^2. Thus its trace is
`-sum_{x != 0} psi(x+a/x^2)`. The map x -> x^-2 permutes the nonzero field
elements. Set h(x^-2)=psi(x) and h(0)=0, then apply the integer Walsh transform.
The output frequency for a has coordinates Tr(a*alpha^j), for 0 <= j < d.
Negating that output gives the trace. No rounding or modular recovery is needed.

Run `make check/native_tables_magma.m` to compare every table entry in both
characteristics through degree 12 with the saved
[prototype data](../tests/native/data/trace_tables/README.md).
The independent check counts all affine pairs (x, y) at every parameter,
including zero, in degrees 1 through 4. It uses finite-field arithmetic without
character or transform tables. Each test run emits a SHA-256 hash covering the
characteristic, degree, polynomial, and complete ordinary table in the fixture
format. These hashes are checked against the retained fixture hashes.

## Native point counts

Build with `make census-native`. The counting command reads one `p degree parameter`
row per line from standard input:

```sh
printf '3 2 0\n3 2 zero\n' | artifacts/native/census-native count
```

An integer parameter e denotes a=alpha^e, with 0 <= e < p^degree-1. The token
`zero` denotes a=0. The CSV output gives `p,degree,parameter,trace,points`, where
points=p^degree+1-trace includes the point at infinity, also on the singular cubic.
The field polynomials, generators, and table order match `FibrationTraceTable`.

`native/field.cpp` retains the prototype's polynomial and packed-coordinate
arithmetic. `native/lfun.cpp` constructs the F_3 tables by exact modular convolution
and the F_2 tables by the integer Walsh transform. Construction checks the prime,
cube root, normalization, coefficient bounds and nodal trace. All modular products
and butterfly sums fit in signed 64-bit integers. Each field and degree is computed
once per process; subsequent counts use its table.

The optional `--table-cache DIRECTORY` writes `f2-d1.txt` through `f3-d12.txt`
as needed, using the canonical three-line fixture format. An existing file must
match the exact regenerated table, including its field polynomial. A changed or
incomplete file causes an error before a count using that table is emitted.
Verification still computes the exact transform; the disk cache records reproducible
tables and does not reduce the arithmetic cost. The command needs no Magma subprocess.

Run `make check/native_tables.m` for fresh Magma comparisons of every native entry
through degree 12 in both characteristics. It also checks direct affine-pair counts
through degree 3, reproducible cache hashes and rejection of changed entries.
