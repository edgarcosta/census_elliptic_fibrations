# Census of elliptic K3 surfaces over F2

This Magma code computes the 1,397,598 isomorphism classes of elliptic K3 surfaces with a section over $\mathbf F_2$ (zero section fixed, automorphisms of the base allowed) in under two minutes.

## What is counted

Each surface has a Weierstrass model $y^2 + a_1xy + a_3y = x^3 + a_2x^2 + a_4x + a_6$ in $\mathbf P(1,1,4,6)$. Here $a_i$ is a binary form of degree $2i$ in $s,t$ over $\mathbf F_2$. Models are identified under $\mathrm{GL}(2,\mathbf F_2)$ acting on $(s,t)$ and the changes $x \mapsto x+q$, $y \mapsto y+px+r$, with $p,q,r$ forms of degree 2, 4, 6. A model with $\Delta \ne 0$ that is minimal at every place is an elliptic K3 surface.

## Method

In characteristic 2 the change $x \mapsto x+a_2$ kills $a_2$, leaving 32 coefficient bits. The rest is the staged quotient and stabilizer enumeration of Lemma 1.1 and Algorithm 1.2 in Auel, Kulkarni, Petok and Weinbaum, "A census of cubic fourfolds over $\mathbf F_2$", [arXiv:2306.09908](https://arxiv.org/abs/2306.09908). It first finds the 88 orbits of pairs $(a_1,a_3)$. For each pair, its stabilizer acts affinely on the quotient of the $(a_4,a_6)$ space by the shifts $y \mapsto y+r$, which has at most $2^{16}$ points. Non-minimality at a place forces $v(\Delta) \ge 12$ and $\deg \Delta = 24$, so minimality only needs checking at $s$, $t$, $s+t$ and $s^2+st+t^2$.

## Files

- `elliptic_fibration_orbits_f2.m`: the orbit enumeration. `FibrationOrbitPrefixesF2` lists the 88 pairs $(a_1,a_3)$, `FibrationFiberOrbitKeysF2` returns one key per orbit over a pair (dropping $\Delta = 0$ by default), and `FibrationFromKeyF2` turns a key back into a model.
- `k3_minimality_f2.m`: `IsEverywhereMinimalF2`, which returns whether a model has $\Delta \ne 0$ and is minimal at the four candidate places, and the places where it is not minimal.
- `census_driver.m`: runs both over all 88 pairs and prints the counts per pair, the census totals and the CPU times.

## Usage

Run `magma -b census_driver.m` from the repository root. It is tested with Magma V2.30-1 and takes about 100 s CPU on one core.

## Results

| Models | Orbits |
| --- | ---: |
| All | 1,402,468 |
| $\Delta = 0$ | 736 |
| $\Delta \ne 0$ | 1,401,732 |
| Non-minimal at exactly one place | 4,129 |
| Non-minimal at two places | 5 |
| Elliptic K3 (everywhere minimal) | 1,397,598 |
| K3 with $a_1 = 0$ ($j = 0$) | 173,977 |
| K3 with $a_1 = s^2$ | 524,493 |
| K3 with $a_1 = st$ | 523,262 |
| K3 with $a_1 = s^2+st+t^2$ | 175,866 |

## Checks

- The 88 pairs, the 1,402,468 orbits and the 736 with $\Delta = 0$ agree with an independent Burnside count over all $2^{37}$ models.
- The minimality filter agrees with Magma's Tate algorithm (`LocalInformation`) on all 97,036 census representatives with $v(\Delta) \ge 12$ at a candidate place, and on 20,000 random others.
- Not yet done: a comparison with an outside count, and the automorphism groups needed for a mass formula check.

Author: Edgar Costa
