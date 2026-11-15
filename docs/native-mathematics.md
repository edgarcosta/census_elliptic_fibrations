# Exact trace tables and the Frobenius-twist sign

This note proves the trace-table recovery bound and the functional-equation
sign formula for the smooth, everywhere minimal elliptic K3 census entries.

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

## Frobenius-twist functional-equation sign

Let K = F_q(t), q = p^f. Let E/K be the generic elliptic curve of a
smooth, relatively minimal elliptic K3 surface X with its zero section.
Here X is the smooth surface obtained by resolving its Weierstrass model.
Put E' = E^(p), with Weierstrass coefficients a_i(t)^p, and let X' be its
smooth relatively minimal elliptic surface. X' need not be K3.
The sign used to complete L(E,T) is

```text
epsilon = (-1)^(f*(h(E')-h(E)) + v_p(C(E')/C(E))),
h(A) = (1/12)*sum_v d_v*v(Delta_min,A),
C(A) = product_v c_v(A),             d_v = [k_v:F_q].
```

Valuations v have v(uniformizer) = 1; v_p is the integer p-adic valuation,
extended to nonzero rational numbers. Each v is a closed place, including
infinity. The Tamagawa number c_v counts k_v-rational points of the Neron
component group. It occurs once in C, without exponent d_v. Discriminant
degrees do carry d_v. In the implemented prime-field case f = 1 and h(E) = 2.

### Minimal models and the height convention

Use the Hodge line bundle of invariant differentials on P^1 for each
minimal elliptic surface. The sections rule out multiple fibres.
Its degree is h(A) = chi(O_XA), and its twelfth
power has the minimal discriminant as a section. Thus its degree gives the
displayed discriminant sum, including wild valuations. These model and
line-bundle conventions are described in [Ulmer's Park City lectures,
Lecture 3, the sections on models and height](https://arxiv.org/html/1101.1939).

For E the census chart has deg a_i <= 2i. At a finite place run Tate
reduction and retain the minimal discriminant valuation and Tamagawa number.
At infinity, with u = 1/t, use x_inf = u^4*x, y_inf = u^6*y and

```text
b_i(u) = u^(2i)*a_i(1/u),             i in {1,2,3,4,6}.
```

The original integral chart has discriminant valuation 24-deg Delta(t).
For E', form a_i(t)^p on the finite chart and b_i(u)^p on the infinity
chart. The latter is u^(2ip)*a_i(1/u)^p: the initial chart weights are 2ip.
Minimize both twisted charts afresh. Translations preserve the discriminant;
each allowed scaling x = pi^2*x_new, y = pi^3*y_new lowers its valuation by
12. Therefore h(E') = p*h(E)-sum_v d_v*s_v, where s_v counts those scalings
from the twisted original minimal models. Setting h(E') = 2p without these
reductions loses the term that can determine the sign.

The relative Frobenius F: E -> E' is a K-isogeny of degree p, with dual
V satisfying VF = [p]. For ell != p it identifies the rational Tate modules,
so L(E,T) = L(E',T) and the conductor exponents agree at every place.
Good reduction is also isogeny invariant. Consequently E and E' have the
same bad-place set, which is sufficient for both discriminant sums and both
Tamagawa products. No new bad places need factoring after the twist.

### BSD normalization and why it is a theorem here

Write r = rank E(K). Use Shioda's rational height pairing, given by negative
intersection of section classes after projection perpendicular to the zero
section and fibre components. Its diagonal convention is

```text
<P,P> = 2*h(A) + 2*(P.O) - sum_v d_v*contr_v(P).
```

The local correction contr_v is per geometric fibre; d_v includes its
conjugates. Let Reg(A) be the determinant on a full basis of A(K)/torsion,
with Reg(A) = 1 in rank zero, and let t_A = #A(K)_tors. The logarithmic
Neron-Tate pairing is log(q) times this pairing. Define

```text
B(A) = lim_(T -> 1/q) L(A,T)/(1-q*T)^r.
```

In this convention the BSD leading-term formula over P^1 is

```text
B(A) = #Sha(A)*Reg(A)*C(A) / (t_A^2*q^(h(A)-1)).       (BSD)
```

To see the power of q, choose a nonzero rational invariant differential
eta. If eta = u_v*eta_v with eta_v a minimal local differential, then
sum_v d_v*v(u_v) = h(A). The differential measure contributes
q*product_v q^(-d_v*v(u_v)) = q^(1-h(A)) over P^1. The remaining local
component factors give C(A). Also 1-q^(1-s) is asymptotic to
log(q)*(s-1), so passing from an s-leading term to B(A) cancels the
factor log(q)^r in the logarithmic regulator. This specializes the
normalization in [Ulmer, Curves and Jacobians over function fields,
"Heights" and "The refined BSD conjecture"](https://arxiv.org/html/1203.5573).
In particular, for the K3 curve (BSD) gives
#Sha(E)*Reg(E) = q*B(E)*t_E^2/C(E), as in the invariant records.

The needed finiteness is unconditional for these entries. The Tate
conjecture holds for X by [Madapusi Pera, Theorem 1](https://arxiv.org/html/1301.6326)
in odd characteristic and [Kim-Madapusi Pera, Theorem 2](https://arxiv.org/html/1512.02540)
in characteristic 2. It gives a finite ell-primary Brauer group for
ell != p. For an elliptic surface with a section, Br(X) is canonically
Sha(E); see [Ulmer's Park City lectures, the Tate-Shafarevich and Brauer
groups and Theorem 40.1](https://arxiv.org/html/1101.1939).
The [Kato-Trihan theorem](https://doi.org/10.1007/s00222-003-0299-2)
then gives finite Sha(E), equality of analytic and algebraic ranks, and
(BSD), including its p-part. Its implication from finiteness of one
ell-primary part is also stated in [Milne's account of the comparison
theorem](https://www.jmilne.org/math/articles/1982f.html).

The maps F and V induce maps on Sha whose composites are multiplication
by p. They therefore identify the ell-primary groups for ell != p.
Thus Sha(E')[ell^infinity] is finite, and Kato-Trihan supplies the full
finiteness and BSD formula for E' as well. This step uses no K3 assertion
about X'. Its zero section gives Br(X') = Sha(E').

### Square order, including the characteristic prime

[Liu-Lorenzini-Raynaud, On the Brauer group of a surface,
Theorem 1, Invent. Math. 159 (2005), 673-676](https://dinolorenzini.franklinresearch.uga.edu/sites/default/files/inline-files/papers/LLR2InventionesPublished.pdf)
proves the required square-order theorem. It applies to a smooth proper
geometrically connected surface over a finite field when one primary part
is finite. There is no restriction excluding p = 2 or the p-primary part.
Both X and X' satisfy the surface hypotheses and have finite Brauer groups
by the preceding argument. Hence #Sha(E) and #Sha(E') are squares.

The [2018 corrigendum, corrected Corollary 3](https://www.math.u-bordeaux.fr/~qliu/articles/CorrigendumToNeronModelsLieAlg.pdf)
retains square order. It adds a square index factor to the general
Brauer-Sha order comparison. Here the zero section makes the global and
local periods and indices all 1. It also identifies Pic^0 of the generic
curve with its Jacobian's rational points. That extra index is therefore 1 too.
The surface theorem supplies the p-part needed here directly; a pairing
argument restricted to primes different from p would not suffice.

### Derivation of the implemented exponent

F is injective on K-points, since its kernel is infinitesimal. On the free
Mordell-Weil groups its image has finite index I, since V composes with F
to [p]. Pullback of the zero divisor is F^*(O') = p*O. Functoriality of
canonical heights, in the same normalization on both curves, gives

```text
<F(P),F(Q)>_E' = p*<P,Q>_E,
Reg(E')/Reg(E) = p^r/I^2.
```

Divide the two instances of (BSD). Their leading terms agree because their
L-functions and ranks agree. With h = h(E), h' = h(E'), C = C(E), C' = C(E'),
the exact identity is

```text
p^r/I^2 = (#Sha(E)/#Sha(E'))*(t_E'/t_E)^2*q^(h'-h)*C/C'.
```

Taking v_p and reducing modulo 2 removes I^2, the torsion square, and both
Sha orders. Therefore

```text
r = f*(h'-h) - v_p(C'/C) mod 2
  = f*(h'-h) + v_p(C'/C) mod 2.
```

Finally, use the sign convention

```text
L(T) = epsilon*(q*T)^D*L(1/(q^2*T)),
c_(D-i) = epsilon*q^(D-2i)*c_i,       c_0 = 1.
```

Write L(T) = (1-q*T)^r*Q(T) with Q(1/q) != 0. Under T -> 1/(q^2*T),
the factor 1-q*T becomes -(1-q*T)/(q*T). Cancel it to order r and let
T -> 1/q. This gives 1 = epsilon*(-1)^r, proving the displayed sign rule.
It is the coefficient-completion sign for L(E,T); it needs no extra
determinant correction for reducible fibre components.

### Connection to the prototype and package local data

In the reference prototype `lfun-csign/tate.cpp`, the lambda `frobenius`
moves the coefficient of t^j to t^(pj). Constants are unchanged because
the input field is F_p. The local Tate loop removes every available scaling.
The `compare` lambda accumulates degree-weighted minimal discriminants and
unweighted v_p differences of Tamagawa numbers. The infinity branch twists
the weight-2i chart before local minimization. Finally, `sum0/12`, `sum1/12`
and `tamagawa_exponent` give exactly h, h' and v_p(C'/C). Its
`twist_height-height+tamagawa_exponent` is the proved exponent for f = 1.

This resolves the normalization and p-primary square-order omissions
identified in `verify-lfun/report.md`. Both files are reference inputs
under the prototype root `artifacts/mw/`, not installed dependencies:

| Reference | SHA-256 |
| --- | --- |
| `lfun-csign/tate.cpp` | `a2bd160ae5b1b5c8847a50361bedbfc8f55adc7c5cd8ae73bec304bebc3c8161` |
| `verify-lfun/report.md` | `2c85b269e41fc3ae091a864f6b63867ce40ce29ca0857030b4dcf4f4ea4d8cea` |

The package's [local records](native-local.md) expose the same inputs through
`FibrationLocalData(E)` and `FibrationLocalData(E : FrobeniusTwist := true)`.
Sum `Degree*DiscriminantValuation` for each height and multiply `Tamagawa`
once per row for C. The existing `make check/native_local_magma.m` controls
include the following; each curve has h(E) = 2:

| Curve | h(E') | C(E) | C(E') | v_p(C'/C) | Exponent | epsilon |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| F_2: y^2+y = x^3+t^6*x | 1 | 4 | 2 | -1 | -2 | +1 |
| F_2: y^2+y = x^3+t^7*x | 4 | 2 | 4 | 1 | 3 | -1 |
| F_3: y^2 = x^3+t^8*x+t+1 | 2 | 1 | 1 | 0 | 0 | +1 |

The test checks reciprocal charts and pins these minimal valuations and
Tamagawa numbers for both curves. In particular, the first twist requires
three local scalings at infinity; using h' = 4 would reverse its
sign. These arithmetic controls check the inputs to the proof. Archived
polynomial agreements and the ten Magma timeouts in the verification report
do not establish the theorem. Weil, functional-equation or Artin-Tate
consistency of a proposed polynomial does not prove its sign.

### Wild conductor and sign are separate computations

The full conductor exponent is f_v = v(Delta_min)+1-m_v, where m_v is the
number of geometric components of the minimal fibre. It equals the tame
exponent (1 for multiplicative, 2 for additive reduction) plus the Swan
contribution. The latter can be nonzero in characteristics 2 and 3.
For these curves D = deg L = sum_v d_v*f_v-4. Discarding wild conductor
terms gives the wrong degree and the wrong amount of point counting.

For example, the first F_2 control has I_9* with 14 components, minimal
discriminant valuation 24 and conductor 11 = 24+1-14, with Swan term 9.
Its twist has III, valuation 12 and conductor 11 = 12+1-2. The F_3 control
has II, valuation 24 and conductor 24, with Swan term 22, on both curves.
The degree-weighted conductor determines D; the height difference and
Tamagawa valuations determine epsilon. Wild reduction is already included
in the minimal local data in this proof. Add no Swan term or separate
local-root-number correction to its exponent.
