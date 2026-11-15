# Exact radix-3 trace tables through degree 12

The characteristic-3 ordinary trace tables have an exact certificate using
convolution modulo a prime M > 2*3^12 with M = 1 mod 3. Every integer coefficient
has absolute value at most 3^d in degree d. This bound gives unique recovery for
all 1 <= d <= 12, including the singular parameter zero.

This proves the table method for stage 1 of [the native plan](plan.md).
It does not implement the package's Magma or native table routines.

## Table and field conventions

The prototype root is `artifacts/mw/`. Its `lfun-counter/field.cpp` and
`lfun-csign/field.cpp` are identical at the inspected snapshot. Each has SHA-256
`43cacd63b6f23bb2402aacaf02c65cfe31c7933ef5ecf3560a990860ed28198c`.
The relevant definitions are `polys3`, `fft3` and the characteristic-3 branch of
`Field::Field`. The encoding operations are in `lfun-counter/field.hpp`.

Fix d with 1 <= d <= 12 and Q = 3^d. Work in

```text
K = F_3[X]/(f_d(X)),     alpha = X mod f_d,
u = u_0 + u_1*alpha + ... + u_(d-1)*alpha^(d-1),
index(u) = u_0 + 3*u_1 + ... + 3^(d-1)*u_(d-1).
```

Each u_i is the integer 0, 1 or 2. A row in `polys3` lists the coefficients
of the monic f_d in ascending powers, including the leading 1:

| d | Coefficients of f_d, constant term first |
| --- | --- |
| 1 | [1, 1] |
| 2 | [2, 2, 1] |
| 3 | [1, 2, 0, 1] |
| 4 | [2, 0, 0, 2, 1] |
| 5 | [1, 2, 0, 0, 0, 1] |
| 6 | [2, 2, 1, 0, 2, 0, 1] |
| 7 | [1, 0, 2, 0, 0, 0, 0, 1] |
| 8 | [2, 2, 2, 0, 1, 2, 0, 0, 1] |
| 9 | [1, 1, 2, 2, 0, 0, 0, 0, 0, 1] |
| 10 | [2, 1, 0, 0, 2, 2, 2, 0, 0, 0, 1] |
| 11 | [1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1] |
| 12 | [2, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1] |

The lower coefficients describe f_d itself. Reduction uses
alpha^d = -sum_(i<d) f_(d,i)*alpha^i. They do not describe the negative of f_d.
Each f_d is irreducible and alpha has multiplicative order Q-1, as checked below.
Thus `exps[e]` is alpha^e, `logs[index(alpha^e)]` is e, and
chi(alpha^e) = (-1)^e for 0 <= e < Q-1. Set chi(0) = 0.

The native bit representation places coefficient 1 in bit i of the low plane,
and coefficient 2 in bit i+16 of the high plane. These planes are disjoint.
`index` converts that representation to the base-3 coordinate index above.
Field addition adds each digit modulo 3, with no carries.
For example, in degree at least 2, index(1+2) is 0, whereas ordinary integer
addition modulo Q gives 1+2 = 3. Convolution therefore uses the additive group
(Z/3Z)^d with this coordinate order. Convolution on the cyclic group Z/QZ
would compute a different table for d > 1.

Magma comparisons must use the specified f_d and map the polynomial residue
class alpha explicitly. The [Magma finite-field creation documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/210)
describes extensions from a supplied irreducible polynomial. A default `GF(3^d)`
may choose another polynomial or generator. Even with `ext<F | f_d>`, the local
degree-1 check returns named generator 1; the root of X+1 must be set to 2.
For d > 1, check f_d(alpha) = 0 and Order(alpha) = Q-1 before comparing arrays.
No compatibility between fields of different degrees is needed for this proof.

## The integer convolution equals the required trace

Define two integer arrays, indexed by the additive coordinates of K:

```text
H(u) = #{x in K : -(x^3+x^2) = u},
J(u) = chi(u),
C(a) = (H * J)(a) = sum_(u in K) H(u)*J(a-u).
```

The prototype builds H by visiting every nonzero x and then adding 1 to H(0)
for x = 0. It builds J(alpha^e) = (-1)^e, leaving J(0) = 0.
In particular, H(0) = 2: both x = 0 and x = -1 map to zero.
Grouping the convolution terms by their preimages gives

```text
C(a) = sum_(x in K) chi(x^3+x^2+a).
```

There are exactly 1+chi(b) solutions of y^2 = b for every b in K.
The projective plane cubic E_a: y^2 = x^3+x^2+a has one point at infinity.
Consequently, its point count, including the singular case a = 0, is

```text
#E_a(K) = 1 + sum_(x in K) (1+chi(x^3+x^2+a))
        = Q+1+C(a).
```

For a != 0 the cubic is nonsingular. The derivative in x is -x, so a singular
affine point would require x = y = a = 0. The point at infinity is smooth.
Thus its elliptic trace is t(a) = Q+1-#E_a(K) = -C(a).
The prototype stores `ordinary[e] = -C(alpha^e)` for 0 <= e < Q-1.
These are logarithmic indices: entry e = 0 means parameter a = 1.
Parameter a = 0 is checked separately and is not an extra filled table entry.

At a = 0 the plane cubic is nodal. Its singular cubic point count yields

```text
C(0) = sum_(x != 0) chi(x^2*(x+1))
     = sum_(x != 0) chi(x+1)
     = sum_(z in K) chi(z) - chi(1)
     = -1.
```

Here sum chi(z) = 0, since multiplication by a nonsquare changes its sign.
Hence -C(0) = 1, agreeing with the prototype's nodal check.
This checks the trace sign at zero. Direct counts at nonzero a also check
the contribution from x = 0.

## Exact transform, cube root and inverse normalization

Let M be prime with M = 1 mod 3 and M > 2*3^12. In R = Z/MZ, choose omega with

```text
omega^3 = 1,     omega != 1.
```

Since R is a field, 1+omega+omega^2 = 0, and omega^-1 = omega^2.
For coordinate vectors k,u in (Z/3Z)^d, use their coordinate dot product
modulo 3 to define

```text
F(v)(k) = sum_u v(u)*omega^(k dot u),
F^-1(w)(a) = Q^-1 * sum_k w(k)*omega^(-k dot a).
```

This dot product refers to the polynomial coordinates, not the field trace
pairing. Both give character parametrizations, but the radix-3 loops use the
coordinate dot product. No trace-dual change of basis is required here.

For each coordinate, sum_(k_i=0..2) omega^(k_i*b_i) is 3 if b_i = 0 and is 0
otherwise. Multiplying these identities gives

```text
sum_k omega^(k dot b) = Q if b = 0, and 0 otherwise.
```

As M != 3, Q is invertible in R. The orthogonality identity proves inversion.
Expanding the product of two forward transforms also proves

```text
F(H * J)(k) = F(H)(k)*F(J)(k),
C mod M = F^-1(F(H mod M)*F(J mod M)).
```

Indeed, writing a = u+v turns omega^(k dot a) into
omega^(k dot u)*omega^(k dot v). This argument uses addition in (Z/3Z)^d.

One coordinate transform takes a triple (x,y,z) to

```text
(x+y+z, x+omega*y+omega^2*z, x+omega^2*y+omega*z).
```

Apply this to triples at positions block+j, block+j+s, block+j+2s.
Use strides s = 1,3,...,3^(d-1), with block a multiple of 3s and 0 <= j < s.
These are the independent coordinate transforms of the tensor product.
There are no additional twiddle factors or digit reversals.
For the inverse, replace omega by omega^2 and multiply every output by Q^-1
once after all d coordinate transforms. Equivalently, normalize each inverse
coordinate transform by 3^-1. Applying both normalizations would be wrong.

The prototype's forward complex root has negative imaginary part and its
inverse root is the conjugate. Either exact root of order 3 can serve as the
forward omega, provided the inverse uses its reciprocal. Swapping the roots
permutes frequencies by k -> -k and leaves the convolution unchanged.

All butterfly and pointwise operations are performed modulo M. Encode -1 as
M-1, use omega^2 in place of complex conjugation, and use the modular inverse
of Q in place of division in floating point. An element of order Q is not
required: every nontrivial character of (Z/3Z)^d has order 3.

## Coefficient bounds and unique integer recovery

H(u) counts preimages, so H(u) >= 0 and sum_u H(u) = Q.
Each H(u) is at most 3 because x^3+x^2+u is a nonzero polynomial of degree 3.
Every J(u) lies in {-1,0,1}. Therefore, for every parameter a,

```text
|C(a)| = |sum_u H(u)*J(a-u)|
       <= sum_u H(u)*|J(a-u)|
       <= sum_u H(u)
       = Q = 3^d.
```

The same bound holds for t(a) = -C(a). This is an elementary bound; it requires
neither a Weil bound nor numerical error estimates.

Reduction modulo M is injective on the integer interval [-Q,Q]. If A and B
in that interval have the same residue, then M divides A-B, while
|A-B| <= 2Q < M. Thus A-B = 0.
For a canonical residue r in [0,M-1], recover the unique centered integer by

```text
center(r) = r       if r <= (M-1)/2,
            r-M    otherwise.
```

The true C(a) is in [-Q,Q] inside that centered interval. Therefore
center(C(a) mod M) = C(a), with no rounding. The recovered trace is its negative.
If a candidate integer trace also satisfies |t_candidate(a)| <= Q, agreement
modulo M proves exact integer equality by the same argument.
Without a bound on the candidate, modular agreement alone is insufficient.

For the existing tables, center the exact convolution, negate it and compare
as signed integers to each widened `ordinary[e]`. Check the recovered integer
against the destination type before narrowing it. Equality of centered
integers certifies the stored value even if the floating construction produced
an incorrect intermediate value. The distance-to-integer tolerance 1e-5 is
not part of this certificate.

One admissible numerical witness is

```text
M = 1062907,          (M-1)/3 = 354302,
omega = 2^354302 mod M = 662797,
omega^2 mod M = 400109,
omega^3 mod M = 1,
1 + 662797 + 400109 = 1062907 = M.
```

Primality is checked by trial division: floor(sqrt(M)) = 1030 and no integer
from 2 through 1030 divides M. Also M mod 3 = 1 and

```text
2*3^12 = 1062882 < 1062907 = M,
3^12 = 531441 < 531453 = (M-1)/2.
```

The following table gives every degree bound and an explicit inverse for the
normalization. In each row M is strictly larger than the displayed 2Q, and
Q times the last entry is 1 modulo M.

| d | Q = 3^d | 2Q | Q^-1 mod 1062907 |
| --- | ---: | ---: | ---: |
| 1 | 3 | 6 | 708605 |
| 2 | 9 | 18 | 590504 |
| 3 | 27 | 54 | 551137 |
| 4 | 81 | 162 | 892317 |
| 5 | 243 | 486 | 297439 |
| 6 | 729 | 1458 | 807751 |
| 7 | 2187 | 4374 | 977855 |
| 8 | 6561 | 13122 | 680254 |
| 9 | 19683 | 39366 | 935356 |
| 10 | 59049 | 118098 | 1020390 |
| 11 | 177147 | 354294 | 340130 |
| 12 | 531441 | 1062882 | 467679 |

The bound is on the final integer convolution, not its transformed factors.
Those factors are residues and need not have small integer lifts.
For this witness, a product of canonical residues is less than
M^2 = 1129771290649 < 2^41, so unsigned 64-bit multiplication suffices.
Reduce products before summing a butterfly. Each sum is then at most
3*(M-1) = 3188718 < 2^22. For another modulus, establish its machine bounds
separately; the mathematical recovery bound remains M > 2*3^12.

## Independent small-field checks and scope of the evidence

The independent oracle enumerates the solutions of
y^2 = x^3+x^2+a directly using Magma's finite-field arithmetic. It adds the
single point at infinity and returns Q+1 minus that count.
It uses neither a Fourier transform nor the prototype's logarithm, quadratic
character or trace tables. Sharing the specified field polynomial and the
coordinate map is necessary to compare the same array entries.

Check every a, including zero, in degrees 1 through 4. These 3+9+27+81 = 120
parameters cover ordinary elliptic traces, the nodal check and both character
signs. They also exercise the least significant digit convention and additions
with no carries. Compare to exact modular convolution and to the existing table at
e whenever a = alpha^e. A direct count also catches a shared sign error that
comparison between two implementations of the same transform would miss.

The scratch mathematical check [trace-proof-check.m](trace-proof-check.m)
was run with `magma -bn`; its first line is `SetQuitOnError(true);` and its
last line is `quit;`. Its [log](trace-proof-check.log) records:

- Primality, root identities, all 12 inequalities and all 12 inverse constants.
- Irreducibility and primitivity of all 12 source polynomials in each of
  characteristics 2 and 3, using the explicit polynomial root in degree 1.
- Exact forward/inverse normalization and direct affine point-count agreement
  at all 120 parameters in degrees 1 through 4, including -C(0) = 1.

The characteristic-2 polynomial checks confirm the adjacent prototype field
convention. Its Walsh transform is integer arithmetic and is outside this
radix-3 proof.

The scratch check does not compare the existing floating tables. The package
Magma computation and native implementation remain separate steps. They must
compare every table entry through degree 12, retaining the field definitions
and reproducible hashes. Direct small-field counts check the construction
independently; the proof above supplies exact recovery for every supported
degree. This document makes no claim that the future implementations or
the existing floating tables have already passed that full comparison.
