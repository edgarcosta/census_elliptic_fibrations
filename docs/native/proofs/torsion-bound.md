# Complete torsion search from the K3 fibre corrections

For every supported smooth, relatively minimal elliptic K3 surface with
section, the sum Cmax of degree-weighted local maximum corrections is at
most 6. Every nonidentity torsion section therefore has (P.O) <= 1.
Here the constant field is F_2 or F_3 and the base is P^1. The generic
fibre must be elliptic: the quasi-elliptic and nonminimal census cells are
outside this assertion. No restriction on the j-map or its separability
is needed.

## Geometric components and their corrections

At a closed place v put d_v = [k_v:F_q] and delta_v = v(Delta_min).
Let m_v count the components after extending k_v to its algebraic
closure. Use one geometric fibre above v when computing its
correction. Its d_v conjugates account for the degree weight below.

A section has intersection 1 with the whole geometric fibre. Consequently,
it meets exactly one component, transversely at a smooth point, and that
component has multiplicity 1. Components of larger multiplicity and
intersections between components cannot be met by a section.
Label the component met by O as 0; its correction is 0.
For any other component i met by P, the correction is (M_v^-1)_(i,i),
where M_v is minus the intersection matrix of the nonidentity components.
These are the positive Cartan matrices of the indicated root lattices.
This is [Shioda's height formula, Theorem 8.6 and table (8.16)](https://rikkyo.repo.nii.ac.jp/record/10027/files/AA00610867_39-02_05.pdf).

The following table includes every Kodaira type over a perfect residue
field, hence every type in characteristics 2 and 3. The column e_v is
m_v for multiplicative reduction and m_v+1 for additive reduction.
For good reduction put e_v = 0.

| Type | Root lattice | m_v | e_v | Corrections on multiplicity-1 components | Geometric maximum c_v^geom |
| --- | --- | ---: | ---: | --- | --- |
| I_0 (good) | none | 1 | 0 | 0 | 0 |
| I_n, n >= 1 | A_(n-1) | n | n | i*(n-i)/n, 0 <= i < n | floor(n^2/4)/n |
| II | none | 1 | 2 | 0 | 0 |
| III | A_1 | 2 | 3 | 0, 1/2 | 1/2 |
| IV | A_2 | 3 | 4 | 0, 2/3 | 2/3 |
| I_n*, n >= 0 | D_(n+4) | n+5 | n+6 | 0, 1, (n+4)/4 | (n+4)/4 |
| IV* | E_6 | 7 | 8 | 0, 4/3 | 4/3 |
| III* | E_7 | 8 | 9 | 0, 3/2 | 3/2 |
| II* | E_8 | 9 | 10 | 0 | 0 |

For I_1 the root lattice is empty. The I_0* nonidentity simple
components all have correction 1. For I_n* with n > 0, the near
simple component has correction 1 and the two far simple components
have correction (n+4)/4. The other vertices of D_(n+4) have fibre
multiplicity 2; their inverse-Cartan diagonal entries are not section
corrections. II* has only the identity component of multiplicity 1.

For example, the A_(n-1) inverse diagonal is i*(n-i)/n.
To check this, solve with boundary values x_0=x_n=0:

```text
2*x_j-x_(j-1)-x_(j+1) = 1 if j=i, and 0 otherwise.
```

The solution has x_i = i*(n-i)/n.
The simple outer vertices of D_(n+4) give 1 and (n+4)/4.
Inverting the E_6 and E_7 Cartan matrices at their nonidentity
multiplicity-1 vertices gives 4/3 and 3/2. Thus wild reduction does
not call for a new height-correction table: the intersection diagrams
and multiplicities still determine these values.

## Wild conductor, discriminant, and infinity

Let w_v be the Swan conductor. Ogg's formula and the tame conductor give

```text
f_v = delta_v + 1 - m_v,
f_v = 0 (good), 1 (multiplicative), 2+w_v (additive),
delta_v = e_v+w_v,       w_v >= 0.
```

Here w_v = 0 for good or multiplicative reduction. It can be positive
at additive fibres in characteristics 2 and 3. These relations and the
unchanged Kodaira classification are given in [Schuett-Shioda,
Sections 4.1, 4.6 and 5.9](https://arxiv.org/html/0907.0298#S4.SS6).
In particular, delta_v is not determined by the Kodaira symbol alone.

Every table row satisfies c_v^geom <= e_v/4. For I_n this follows
from i*(n-i) <= n^2/4. For I_n* it is (n+4)/4 <= (n+6)/4.
The remaining positive maxima satisfy 1/2 <= 3/4, 2/3 <= 4/4,
4/3 <= 8/4, and 3/2 <= 9/4; the zero maxima satisfy the inequality too.
Consequently

```text
c_v^geom <= e_v/4 <= delta_v/4.
```

The extra wild valuation only increases the discriminant cost of a given
component correction. Replacing delta_v by the tame e_v in the global
discriminant sum would discard that cost.

The Hodge bundle has degree chi(O_X) = 2, and its twelfth power has
the minimal discriminant as a section. Thus

```text
sum_v d_v*delta_v = 12*chi(O_X) = 24.
```

Equivalently, Ogg's formula gives
deg(conductor) + sum_v d_v*(m_v-1) = 24, including all wild terms.

See [Ulmer, Section 34, the Hodge bundle and height](https://arxiv.org/html/1101.1939#S34)
for the bundle and discriminant convention. The sum runs over all closed
places, including infinity. In the census chart, infinity has degree 1
and coefficients u^(2i)*a_i(1/u), with u = 1/t. Its minimal discriminant
valuation is 24-deg Delta(t) for an everywhere minimal K3 equation.
Apply the same Kodaira and wild-conductor formulas to this local chart.
There is no separate bound or omitted correction at infinity.

## Residue-field rational components

A section over F_q(t) specializes to a k_v-rational point. Its geometric
component must therefore be fixed by Frobenius of k_v. Here Frobenius is
the q^d_v power map on one geometric fibre; q-Frobenius also permutes
the d_v conjugate fibres. Only multiplicity-1 fixed components can
contribute. The following restrictions follow from the automorphisms of
the Kodaira diagrams fixing the identity component:

- For split I_n all indices i are allowed. For nonsplit I_n Frobenius
  sends i to -i modulo n. Only i=0 and, for even n, i=n/2 are fixed.
  The nonsplit maximum is therefore 0 for odd n and n/4 for even n.
- For I_n* with n > 0, the near simple component is fixed. The far pair
  is either fixed individually (Tamagawa number 4) or exchanged
  (Tamagawa number 2). The respective maxima are (n+4)/4 and 1.
  For I_0*, Frobenius permutes the three nonidentity simple components.
  A 3-cycle leaves only the identity available (Tamagawa number 1);
  a transposition or the identity gives maximum 1 (number 2 or 4).
- For IV and IV*, the two nonidentity simple components are both fixed
  or exchanged. Tamagawa number 3 allows the displayed positive
  correction; number 1 allows only correction 0.
- For III and III*, the unique nonidentity simple component is fixed.
  For II and II*, only the identity simple component is available.

These are necessary component restrictions; they do not assert that every
choice occurs for a global section. In particular, do not average a
correction over an orbit of components or use an orbit sum as a component
that a section can meet.

Let c_v^rat be the largest correction among these allowed components,
including the identity, and define

```text
Cmax_rat = sum_v d_v*c_v^rat,
Cmax_geom = sum_v d_v*c_v^geom.
```

These are sums of local maxima; a section need not attain them
simultaneously. Either convention gives a valid upper bound Cmax:

```text
sum_v d_v*contr_v(P) <= Cmax_rat <= Cmax_geom
                   <= (1/4)*sum_v d_v*delta_v = 6.
```

## The complete section bound and the identity

Use the height normalization from the [sign proof](frobenius-sign.md). For P != O,
the intersection m = (P.O) is a nonnegative integer and

```text
<P,P> = 4 + 2*m - sum_v d_v*contr_v(P).
```

The height pairing has kernel exactly the torsion subgroup, including
characteristic-primary torsion. This follows from its intersection-theoretic
construction; see [Shioda, Theorem 8.4](https://rikkyo.repo.nii.ac.jp/record/10027/files/AA00610867_39-02_05.pdf).
It also follows over the finite constant field from the canonical height
and finiteness of points of bounded height. Thus, every nonidentity
zero-height section satisfies

```text
4+2*m = sum_v d_v*contr_v(P) <= Cmax <= 6,
0 <= m <= floor((Cmax-4)/2) <= 1.
```

If Cmax < 4 there are no nonidentity torsion sections. Otherwise, an
exhaustive search through that integer bound, followed by an exact
zero-height check, recovers every nonidentity torsion point. This argument
does not use specialization injectivity on p-torsion or a separability
assumption on multiplication by its order.

The identity must be included separately. Its height is 0 but
(O.O) = -chi(O_X) = -2. Substitution in the same height formula gives
0 = 4+2*(-2)-0, while the nonnegative-intersection argument concerns
two distinct sections. O is also absent from affine coordinate searches.

For completeness in the polynomial chart, write a nonidentity point as
x=A/g^2, y=B/g^3 with g monic and gcd(A,g)=1. The formal parameter at
O gives poles of orders 2*r and 3*r at an intersection of multiplicity r.
Thus deg g is the finite part of (P.O). Using x_inf=u^4*x and
y_inf=u^6*y at infinity gives

```text
deg g <= m,       deg A <= 4+2*m,       deg B <= 6+3*m.
```

Enumerating all such triples through m <= 1 includes denominators of
degree 0 or 1 and intersections with O at infinity. Checking the equation,
the actual intersection and the exact height removes extra candidates.
No search beyond m=1 is needed for torsion on the supported K3 population.

The prototype `c-stages/reference.m` sums geometric maxima as
`f[2]*f[6]`, with f[2] the closed-place degree. The function `fibre_data`
in `hard/code/lattice.py` restricts to rational components and returns
corrections already multiplied by the place degree. Its output must not
be weighted a second time. These paths under the prototype root
`artifacts/mw/` fix conventions only. The universal bound follows from
the local inequality and the K3 discriminant degree, independently of
any sampled agreement or prototype torsion computation.
