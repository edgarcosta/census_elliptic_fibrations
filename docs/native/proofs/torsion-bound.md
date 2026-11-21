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

The [Hodge bundle](https://arxiv.org/html/1101.1939#S34) has degree chi(O_X) = 2, and its twelfth power has
the minimal discriminant as a section. Thus

```text
sum_v d_v*delta_v = 12*chi(O_X) = 24.
```

Equivalently, Ogg's formula gives
deg(conductor) + sum_v d_v*(m_v-1) = 24, including all wild terms.

The sum runs over all closed
places, including infinity. In the census chart, infinity has degree 1
and coefficients u^(2i)*a_i(1/u), with u = 1/t. Its minimal discriminant
valuation is 24-deg Delta(t) for an everywhere minimal K3 equation.
Apply the same Kodaira and wild-conductor formulas to this local chart.
There is no separate bound or omitted correction at infinity.

The degree-weighted geometric maximum is
Cmax = sum_v d_v*c_v^geom <= sum_v d_v*delta_v/4 = 6.

## The complete section bound and the identity

For P != O,
the intersection m = (P.O) is a nonnegative integer and the [normalized height](frobenius-sign.md) is

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
