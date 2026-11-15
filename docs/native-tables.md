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
