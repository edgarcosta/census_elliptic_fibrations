# Frobenius-twist functional-equation sign

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

## Minimal models and the height convention

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

## BSD normalization and why it is a theorem here

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

## Square order, including the characteristic prime

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

## Derivation of the implemented exponent

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

## Connection to the prototype and package local data

In the reference prototype `lfun-csign/tate.cpp`, the lambda `frobenius`
moves the coefficient of t^j to t^(pj). Constants are unchanged because
the input field is F_p. The local Tate loop removes every available scaling.
The `compare` lambda accumulates degree-weighted minimal discriminants and
unweighted v_p differences of Tamagawa numbers. The infinity branch twists
the weight-2i chart before local minimization. Finally, `sum0/12`, `sum1/12`
and `tamagawa_exponent` give exactly h, h' and v_p(C'/C). Its
`twist_height-height+tamagawa_exponent` is the proved exponent for f = 1.

The package's [local records](../local.md) expose the same inputs through
`FibrationLocalData(E)` and `FibrationLocalData(E : FrobeniusTwist := true)`.
Sum `Degree*DiscriminantValuation` for each height and multiply `Tamagawa`
once per row for C. The `tests/native_local_magma.m` controls
include the following; each curve has h(E) = 2:

| Curve | h(E') | C(E) | C(E') | v_p(C'/C) | Exponent | epsilon |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| F_2: y^2+y = x^3+t^6*x | 1 | 4 | 2 | -1 | -2 | +1 |
| F_2: y^2+y = x^3+t^7*x | 4 | 2 | 4 | 1 | 3 | -1 |
| F_3: y^2 = x^3+t^8*x+t+1 | 2 | 1 | 1 | 0 | 0 | +1 |

The test checks reciprocal charts and pins these minimal valuations and
Tamagawa numbers for all three curves. In particular, the first twist requires
three local scalings at infinity; using h' = 4 would reverse its
sign. These arithmetic controls check the inputs to the proof.
Weil, functional-equation or Artin-Tate consistency of a proposed polynomial
does not prove its sign.

## Wild conductor and sign are separate computations

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
