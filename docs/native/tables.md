# Exact ordinary trace tables

`FibrationTraceTable(p,d)` gives exact traces for p in [2,3] and 1 <= d <= 12.
The record contains `Polynomial` (the Conway polynomial), `Field`,
`Generator` (a primitive Conway root alpha), `Ordinary`, and `NodalTrace`.
`Ordinary[e+1]` is the trace at alpha^e; `NodalTrace` is 1 at parameter zero.
Over F_3 in degree 1, alpha=2, although Magma's defining polynomial differs.

The curves are y^2+x*y=x^3+a over F_2 and y^2=x^3+x^2+a over F_3.
At zero they are singular projective cubics with one point at infinity.
Tables use O(d*p^d) arithmetic operations and O(p^d) storage.

## Characteristic 2

Write psi(z)=(-1)^Tr(z). For x != 0, division by x^2 gives
(y/x)^2+y/x=x+a/x^2, so the trace is -sum_(x != 0) psi(x+a/x^2).
At x=0 there is exactly one y, contributing zero to the trace.
The map x -> x^-2 permutes the nonzero elements. Set h(x^-2)=psi(x)
and h(0)=0. The integer Walsh transform of h at the coordinates
Tr(a*alpha^j), 0 <= j < d, is the negative trace.
No modular recovery or rounding is needed.

## Characteristic 3

Let Q=3^d, chi(0)=0 and chi(alpha^e)=(-1)^e. On the additive group
of polynomial coordinates (Z/3Z)^d, set

```text
H(u) = #{x in K : -(x^3+x^2) = u},
J(u) = chi(u),
C(a) = (H*J)(a) = sum_(x in K) chi(x^3+x^2+a).
```

There are 1+chi(b) solutions of y^2=b, hence #E_a(K)=Q+1+C(a)
and t(a)=-C(a). The derivative in x is -x, so nonzero a gives a smooth
curve. At a=0, C(0)=sum_(x != 0) chi(x+1)=-chi(1)=-1,
which proves the nodal trace is 1.

Take M=1062907 and omega=662797 in Z/MZ. The integer M is prime,
omega^3=1, omega!=1 and omega^2=400109=omega^-1. Define

```text
F(v)(k) = sum_u v(u)*omega^(k dot u),
F^-1(w)(a) = Q^-1 * sum_k w(k)*omega^(-k dot a).
```

The dot product uses polynomial coordinates. Character orthogonality gives
sum_k omega^(k dot b)=Q when b=0 and zero otherwise. Thus F^-1 is the
inverse, and expanding the sums gives F(H*J)=F(H)*F(J).
Each coordinate transform sends (x,y,z) to
(x+y+z, x+omega*y+omega^2*z, x+omega^2*y+omega*z).
The inverse uses omega^2 and multiplies by Q^-1 once after all coordinates.

Since H>=0, sum H=Q and |J|<=1, every integer C(a) satisfies |C(a)|<=Q.
The bound M>2*3^12=1062882 makes reduction injective on [-Q,Q].
Centering a residue r by subtracting M when r>(M-1)/2 therefore recovers
C(a) exactly; its negative is the trace. Computed bounds and normalization
are checked during construction.

A product of canonical residues is less than M^2<2^41; a butterfly sum
is less than 2*M^2<2^42. The native operations fit in signed 64-bit integers.

## Native point counts and checks

`census-native count` reads `p degree parameter` rows. An integer parameter e
means alpha^e, with 0 <= e < p^degree-1; `zero` means parameter zero.
The CSV columns are `p,degree,parameter,trace,points`, with
points=p^degree+1-trace, including the point at infinity.
Tables are reused within each process.

Both implementations compare every table with the
[canonical digests](../../tests/fixtures/trace_tables/README.md).
Magma's `TraceOfFrobenius` independently checks every nonzero parameter
through degree 4. The recovery bound applies through degree 12.
