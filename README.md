# Census of elliptic K3 surfaces over F2

This Magma code computes the 1,397,598 isomorphism classes of elliptic K3 surfaces with a section over $\mathbf F_2$ (zero section fixed, automorphisms of the base allowed) in under two minutes. The code for $\mathbf F_3$, described under [Over F3](#over-f3), finds 1,307,723,786 such classes there.

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

## Over F3

### Method

Over $\mathbf F_3$, completing the square gives $y^2 = x^3 + Cx^2 + Ax + B$ with $C$, $A$, $B$ binary forms of degree 4, 8, 12. Models are identified under $\mathrm{GL}(2,\mathbf F_3)$ acting on $(s,t)$ and $x \mapsto x+q$ with $q$ of degree 4; unlike over $\mathbf F_2$, $C$ cannot be removed. Scalars act trivially, so the group acting has order $24 \cdot 3^5 = 5832$, and $\Delta = C^2A^2 - A^3 - C^3B$.

The first stage lists the 1,860 orbits of pairs $(C,A)$. There are 933 with $C \ne 0$, where $A$ is taken modulo $C \cdot V_4$ ($V_d$ being the forms of degree $d$), and 927 with $C = 0$. The stabilizer of a pair acts affinely on what is left. When $C \ne 0$ that is all of $V_{12}$, with $3^{13}$ points; when $C = 0$ it is $V_{12}$ modulo the image of $q \mapsto q^3 + Aq$, with $3^8$ or $3^9$ points.

A model is non-minimal at a place $\pi$ exactly when $\pi^2 \mid C$ and some $r$ modulo $\pi^2$ gives $\pi^4 \mid A - Cr$ and $\pi^6 \mid B + Ar + Cr^2 + r^3$. As over $\mathbf F_2$, only places of degree at most 2 can fail, and over $\mathbf F_3$ there are 7 of them. For 1,758 of the pairs no place satisfies the first two conditions, so every model over them with $\Delta \ne 0$ is minimal.

### Files

- `ef_common.m`: binary forms and base-$p$ keys.
- `ef3_models.m`: models $[C,A,B]$, their keys, and the group action.
- `ef3_orbits.m`: `FibrationOrbitPrefixesF3` lists the 1,860 pairs, `EF3_FiberData` gives the affine maps their stabilizer induces, and `FibrationFiberOrbitKeysF3` returns one key per orbit over a pair.
- `k3_minimality_f3.m`: `IsEverywhereMinimalF3`, and `K3F3_PrefixAutomaticallyMinimal` for the pairs over which every smooth model is minimal.
- `ef3_stream.m`: `CensusF3` runs the chosen pairs, writes one CSV row of counts per pair, and writes representatives as 8-byte little-endian keys.
- `ef3_burnside.m`: a Burnside count over all 5832 group elements, independent of the other files.
- `bench_f3.m`: times one complete pair of each kind.
- `checks/`: the check scripts.

### Usage

Load `ef_common.m`, `ef3_models.m`, `ef3_orbits.m`, `k3_minimality_f3.m` and `ef3_stream.m` in that order, then call for example `CensusF3(R, [], "runs/full")` with `R<s,t> := PolynomialRing(GF(3),2)`; an empty list of indices means every pair. Run the checks from the repository root, as in `magma -b checks/check_orbits.m`. A failed assertion prints `Runtime error in assert`, but `magma -b` still exits with status 0, so read the output. `checks/check_stream.m` writes `stream_check_*` files in the current directory and refuses to overwrite them, so delete them before running it again.

### Results

| Models | Orbits |
| --- | ---: |
| All | 1,307,814,237 |
| $\Delta = 0$ | 338 |
| $\Delta \ne 0$ | 1,307,813,899 |
| Non-minimal at exactly one place | 90,105 |
| Non-minimal at two places or more | 8 |
| Elliptic K3 (everywhere minimal) | 1,307,723,786 |
| K3 with $C = 0$ | 5,451,380 |
| K3 with $C \ne 0$ | 1,302,272,406 |

### Checks

These all pass with Magma V2.30-1.

- A Burnside count over all 5832 elements gives 1,307,814,237, and so does the sum of the Burnside counts of the 1,860 stabilizers on their fibres.
- The 338 orbits with $\Delta = 0$ come from a separate argument: 23 with $C \ne 0$, each equivalent to some $(C,0,0)$, and 315 with $C = A = 0$.
- Explicit enumeration agrees with Burnside on 35 complete fibres, one of each kind (C zero or not, fibre dimension, number of distinct maps, and whether every smooth model is minimal).
- The minimality test agrees with Magma's Tate algorithm (`LocalInformation`) on 994 local checks over 142 models.
- The census without representatives, `CensusF3(R, [], stem : WriteRepresentatives := false)`, took 6,567 s of CPU in total, run as 12 processes on disjoint sets of pairs. Writing the 1.3 billion representatives (10.5 GB) adds an estimated 50 minutes.

Author: Edgar Costa
