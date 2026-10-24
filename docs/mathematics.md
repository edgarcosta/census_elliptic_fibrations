# Mathematics of the census

The census counts degree-bounded Weierstrass equations up to changes of
coordinates. It records the generic curve and the local minimality of the
equation as separate properties. The regular, everywhere minimal cases give
elliptic or quasi-elliptic K3 fibrations with a section.

## Equations and equivalence

Let k = F_p, with p = 2 or 3, let K = k(P1), and put
V_d = H^0(P1,O(d)). The equations are

```text
y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6,
a_i in V_(2i).
```

They lie in weighted projective space P(1,1,4,6), with coordinates s,t,x,y.
The equivalence relation preserves the zero section and the coefficient
line bundle O(2). It allows GL(2,k) on s,t and the changes

```text
x_old = u^2*x + r,
y_old = u^3*y + q*u^2*x + v,
u in k*, r in V4, q in V2, v in V6.
```

The implementation applies these formulas to singular equations as well.
In particular, it counts an equation before any local reduction; division
by a local parameter is not part of this equivalence relation.

Over F2, translation in x removes a2. Over F3, completing the square
removes a1 and a3. The normal forms are

```text
F2: y^2 + a1*x*y + a3*y = x^3 + a4*x + a6,
    degrees (a1,a3,a4,a6) = (2,6,8,12).
F3: y^2 = x^3 + C*x^2 + A*x + B,
    degrees (C,A,B) = (4,8,12).
```

Unit scalings add no further identifications on these normal forms over
the two supported fields. Over F3, scalar base matrices act trivially on
the even-degree forms. The [key guide](keys.md#normalization) gives the
normalization formulas and the encoding of the resulting coefficients.

## Prefixes and affine G-modules

The staged orbit calculation first chooses orbits of coefficient pairs,
then counts the orbits in each fibre under its stabilizer. This is the
quotient-and-stabilizer method of [Auel, Kulkarni, Petok and Weinbaum,
Lemma 1.1 and Algorithm 1.2](https://arxiv.org/abs/2306.09908).
The linear actions on binary forms are symmetric powers of the natural
GL(2,k)-module. Residual translations make the fibre actions affine.

Over F2, the prefixes are the orbits of (a1,a3). A quadratic y-shift q
preserves a2 = 0 when accompanied by the quartic x-shift
r = q*a1+q^2. The stabilizer therefore includes the quadratic shifts
that fix the chosen pair, as well as the appropriate base changes.
The sextic shifts have image

```text
V6 -> V8 x V12,
v |-> (a1*v, v^2+a3*v).
```

This is linear over F2. Quotienting by its image leaves a vector space
with an affine stabilizer action. The package has 88 prefixes over F2.

Over F3, the prefixes are (C,A). When C is nonzero, first choose a
base orbit of C, then an orbit in V8/(C*V4) under its stabilizer.
The x-shift needed to preserve A is determined by the base change,
and its action on V12, the space of B, is affine.
When C = 0, choose a base orbit of A and use the quotient

```text
V12 / {r^3+A*r : r in V4}.
```

The map defining this subspace is linear over F3. The stabilizer of A
acts on the quotient. There are 1,860 prefixes over F3.

## Burnside counting and the six cells

For a finite invariant set T in a fibre with affine group E, write

```text
n(T) = #(T/E) = (1/#E) * sum_{g in E} #Fix_T(g).
```

Fixed points on an affine subspace are found by solving linear equations.
In the fibre W, let S be the discriminant-zero locus, Z the nonregular
locus, and N the nonminimal locus. Each is a union of affine subspaces;
Z is contained in S. The regular cuspidal locus is Q = S minus Z.
The three loci and their intersections are E-invariant.

The engine counts unions as sets, either by deduplicating their points
or by signed intersections of affine subspaces. Adding the counts of
individual bad-place subspaces would count their intersections repeatedly.
Burnside's lemma applies to the invariant union; individual components
need not be invariant. Exact membership in S, Z and N gives the six cells:

| Generic-fibre kind | Everywhere minimal | Nonminimal somewhere |
| --- | --- | --- |
| smooth | n(W) - n(S) - n(N minus S) | n(N minus S) |
| regular cusp | n(Q) - n(N intersect Q) | n(N intersect Q) |
| nonregular | n(Z) - n(N intersect Z) | n(N intersect Z) |

Counting does not enumerate all representatives. When representatives are
requested, the engine enumerates orbits and selects them by membership
in these loci, including the wholly singular prefixes.

## Generic-fibre kind

For a normalized F2 equation, the discriminant is

```text
Delta = a3^4 + a1^3*a3^3 + a1^4*a4^2
        + a1^5*a3*a4 + a1^6*a6.
```

For a normalized F3 equation, it is

```text
Delta = C^2*A^2 - A^3 - C^3*B.
```

In either characteristic, Delta != 0 means the generic curve is smooth.
When Delta = 0, regularity over K distinguishes quasi-elliptic curves
from nonregular curves. Regularity over K does not imply geometric
smoothness over K, since K is imperfect.

Over F2, a regular cusp is characterized by

```text
a1 = a3 = 0, and (a4 is not a square in K or a6 is not a square in K).
```

See [Cossec, Dolgachev and Liedtke, Section 4.4, Case 4, pp. 400-401](https://sites.lsa.umich.edu/idolga/wp-content/uploads/sites/1334/2024/08/EnriquesOne.pdf).
In an affine chart, write a4 = f(z), a6 = g(z). At the geometric cusp,
the derivative in z is f'(z)*sqrt(f(z))+g'(z).
It vanishes exactly when f and g are both squares. If f' is nonzero,
the vanishing of f'(z)*sqrt(f(z))+g'(z) would put sqrt(f) in K,
contradicting f' != 0.
If f' = 0, the condition is g' = 0.
If a1 != 0 and Delta = 0, the singular point is K-rational and nodal.
The case a1 = 0, a3 != 0 is smooth because Delta = a3^4.

Over F3, a regular cusp is characterized by

```text
C = A = 0, and B is not a cube in K.
```

See [Lang, Section 1A, Proposition 1.3, p. 478](https://www.numdam.org/article/ASENS_1979_4_12_4_473_0.pdf).
Here B is a cube exactly when its derivative in an affine chart vanishes.
If C != 0 and Delta = 0, the generic curve has the K-rational node
(x,y) = (A/C,0). If C = 0 and A != 0, Delta = -A^3 != 0.

Squares and cubes here are taken in K, before algebraic closure.
For these forms over a perfect constant field, they are equivalent to
homogeneous polynomial powers. The remaining discriminant-zero equations
are nonregular: nodal curves or cusps with K-rational cusp.

## Local minimality

Represent a closed place by a homogeneous irreducible form pi.
An equation is nonminimal there if an integral translation makes pi^i
divide a_i for i = 1,2,3,4,6. Then x = pi^2*X and y = pi^3*Y lower the
coefficient line bundle. This is the weighted reduction criterion in
[Schuett and Shioda, Section 4.8](https://arxiv.org/html/0907.0298v3).
The same divisibility property is used for nonregular equations; it
does not assert that their resolutions are genus-one fibrations.

Over F2, after normalization, put x_old = x+r and y_old = y+q*x+v.
The transformed coefficients are

```text
b1 = a1,
b2 = r + q*a1 + q^2,
b3 = a3 + r*a1,
b4 = a4 + q*a3 + (v+r*q)*a1 + r^2,
b6 = a6 + r*a4 + r^3 + v*a3 + r*v*a1 + v^2.
```

It suffices to take q modulo pi, r modulo pi^2 and v modulo pi^3.
Choose r = q*a1+q^2 to make b2 zero. For fixed q,v, the divisibility
conditions on b4,b6 cut out an affine subspace in V8 x V12.
See also the quasi-elliptic local reduction in
[Cossec, Dolgachev and Liedtke, (4.4.37)-(4.4.38), pp. 415-416](https://sites.lsa.umich.edu/idolga/wp-content/uploads/sites/1334/2024/08/EnriquesOne.pdf).

Over F3, completing the square is integral. It remains to find a quartic
r modulo pi^2 satisfying

```text
pi^2 | C,
pi^4 | A-C*r,
pi^6 | B+A*r+C*r^2+r^3.
```

For fixed r the last condition describes a coset of
pi^6*V_(12-6*deg(pi)). For the regular cuspidal case, compare
[Lang, Section 1B, pp. 479-480](https://www.numdam.org/article/ASENS_1979_4_12_4_473_0.pdf).

These tests use coefficient divisibility, including when Delta = 0.
A discriminant valuation bound alone would not test those equations.

### Candidate places, including infinity

Outside the cuspidal prefix, the first nonzero coefficient bounds the
degree of a bad place. Over F2, pi divides a1, or pi^3 divides a3
when a1 = 0. Over F3, pi^2 divides C, or pi^4 divides A when C = 0.
Their degrees force deg(pi) <= 2.

In the cuspidal prefixes, rational places suffice unless every place is
bad. Here is the bound without using a nonzero discriminant.
Over F2, r = q^2. The differentials of a4 and a6 are sections of
O(6) and O(10). If da4 != 0, a bad place requires pi^4 to divide da4,
so 4*deg(pi) <= 6. If da4 = 0 and da6 != 0, it requires pi^6 to
divide da6, so 6*deg(pi) <= 10. Both force a rational place.

If a4 = f^2 and a6 = g^2, with f of degree 4, a bad place requires
pi^2 to divide f+q^2. When df != 0, it is a section of O(2), again
forcing deg(pi) = 1. When df = 0, a4 is a fourth power. Global x
and y translations remove a4 and a6, giving y^2 = x^3.

Over F3, for a regular cusp dB != 0 is a section of O(10).
The bad-place equation B+r^3 = pi^6*h gives dB = pi^6*dh,
so 6*deg(pi) <= 10. If B is a cube, translate it away globally.

Thus, the finite candidate lists, shown without specifying their order, are

```text
F2: s, t, s+t, s^2+s*t+t^2.
F3: s, t, s+t, s-t, s^2+t^2, s^2+s*t-t^2, s^2-s*t-t^2.
Cuspidal prefixes: only the linear forms in these lists.
```

The homogeneous conditions include t = 0, the place at infinity in the
chart a_i(s,1). The all-places case is recognized algebraically first:
over F2 it has a1 = a3 = 0, a4 a fourth power and a6 a square.
Over F3 it has C = A = 0 and B a cube. These equations are globally equivalent
to y^2 = x^3. Other equations have a finite bad-place list, possibly empty.
The [interface](interface.md#properties-and-minimality) specifies the certificates.

## Why the regular minimal cases are K3

Resolving a minimal Weierstrass equation with regular generic curve gives
a genus-one fibration with a section and fundamental bundle O(2).
The section rules out multiple fibres. The canonical bundle formula gives

```text
omega_X = f^*(omega_P1 tensor O(2)) = O_X,
R^1 f_* O_X = O(-2),
H^1(X,O_X) = 0.
```

These are the K3 conditions, for both the smooth and regular cuspidal
generic curves. See [Cossec, Dolgachev and Liedtke, Proposition 4.3.7,
p. 385](https://sites.lsa.umich.edu/idolga/wp-content/uploads/sites/1334/2024/08/EnriquesOne.pdf)
for the genus-one canonical bundle formula, and
[Schuett and Shioda, Sections 2.5, 4.8 and 4.10](https://arxiv.org/html/0907.0298v3)
for the elliptic Weierstrass construction and local reduction.

## Counts and independent checks

The [README table](../README.md#results) records the six cells.
Running NumberOfFibrations over each field gives:

| Quantity | F2 | F3 |
| --- | ---: | ---: |
| All orbits | 1,402,468 | 1,307,814,237 |
| Smooth, everywhere minimal | 1,397,598 | 1,307,723,786 |
| Smooth, nonminimal somewhere | 4,134 | 90,113 |
| Regular cusp, everywhere minimal | 680 | 294 |
| Regular cusp, nonminimal somewhere | 38 | 20 |
| Nonregular, everywhere minimal | 7 | 12 |
| Nonregular, nonminimal somewhere | 11 | 12 |

The fast suite checks these counts and explicit equations through the public
interface. The [reference suite](../reference/README.md) independently
enumerates the singular orbits, tests their local congruences in both
charts, and compares the package with scan and Burnside calculations.
The F2 smooth comparison covers every prefix. Over F3, every prefix total
is checked against an independent Burnside count. The smooth count is
checked prefix by prefix on the 930 scanned prefixes, and otherwise only
summed over all prefixes. The smooth minimal/nonminimal split is sampled.
Tate's algorithm supplies further smooth local checks.
