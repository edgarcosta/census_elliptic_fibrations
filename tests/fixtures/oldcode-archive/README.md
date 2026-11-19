# Archived coordinate-orbit multisets

These fixtures compare the package with a completed run of the old census.
Thanks for helping preserve an independent check of the census data.
Only per-prefix counts and mathematical digests are stored here, in Git LFS.
The sorted orbit lists stay on the compute host for diagnosing mismatches.

## Provenance and selection

The old source is [census_elliptic_fibrations at the pinned revision](https://github.com/edgarcosta/census_elliptic_fibrations/tree/2f383c61d15c4ce8b36bc4c8b06d96a9fc173074).
The archive was completed and validated on the EPYC on 2026-11-15, using
Magma 2.29-10. Its root there is
`/home/edgarcosta/projects/census_elliptic_fibrations/_worktrees/epyc-oldcode-2f383c6-20261115/_worktrees/old/_worktrees/old-archive`.
Both launch and validation exited zero, checking all 1,948 prefixes against
the pinned totals: 1,402,468 F_2 orbits and 1,307,814,237 F_3 orbits.

The decoder-only extraction completed on the same day in 207 wall seconds.
Its output is under
`/home/edgarcosta/projects/census_elliptic_fibrations/_worktrees/archive-tests-v1-20261115/results`
on that host, with `fixture/summary.json` and `fixture/orbits.bin` in each
`f2-pNNNN` or `f3-pNNNN` directory. Extraction decoded saved keys and applied
the independent canonicalization below; it did not rerun the old census.
The old source stays outside this repository.

`summaries.tsv` preserves the transferred counts and digests, with these
complete-prefix selections:

| Field | Old prefix indices | Orbits | Smooth minimal | Smooth nonminimal | Delta=0 |
| --- | --- | ---: | ---: | ---: | ---: |
| F_2 | 1-88 | 1,402,468 | 1,397,598 | 4,134 | 736 |
| F_3 | 1, 934 | 1,594,638 | 1,594,323 | 0 | 315 |

The selected F_3 prefixes have canonical labels 0 and 21870, represented by
(C,A) = (0,0) and (t^4,s^7*t). This bounded selection exercises C=0 and
C!=0; it does not sample all F_3 stabilizers or smooth-nonminimal models.
The archive's false minimality flag at Delta=0 does not classify singular
minimality. Its four singular cells are compared only through their sum.
The separate classification, Tate, and formula tests check finer cells.

The transferred file has SHA-256
`54c3ec33192f354d3f9f84a456d3b36a9d9922f98140ec7d2efd6aa29499afdf`.
This records transfer integrity only; tests never use a file-byte checksum
as a mathematical oracle. There are 360 prefix/bucket rows, followed by
`CELLS`, eight combined rows, and `END_SUMMARIES`. The old ordinal is metadata;
comparisons identify prefixes by their canonical mathematical labels.

## Complete invariant and canonicalization

For each equation, I denotes the lexicographically least normalized coefficient
tuple in its allowed degree-bounded coordinate orbit. Coefficients of a
degree-d form are ordered s^d, s^(d-1)t, ..., t^d, including zero coefficients.
Over F_2 the tuple concatenates a1,a3,a4,a6, with lengths 3,7,9,13.
Over F_3 it concatenates C,A,B, with lengths 5,9,13.
An integer packs this tuple in base q, most significant digit first.
This integer encodes the canonical form independently of the input representative key.
The prefix label packs its first 10 or 14 digits in the same way.

Over F_2, enumerate GL(2,2) base changes and degree-2 forms h.
After the base change, use x_old=x+r, y_old=y+h*x+z, where
r=h^2+a1*h and z has degree 6. The condition on r preserves a2=0.
Retain every transport attaining the least prefix (a1,a3).
At that target, the z action adds (a1*z,z^2+a3*z) to (a4,a6).
Reduce suffixes modulo this linear subspace and take their least value.

Over F_3, use the 24 projective base matrices whose first nonzero entry is 1.
Scalar base matrices and Weierstrass unit changes act trivially on these forms.
A quartic translation h sends
(C,A,B) to (C,A-C*h,B+A*h+C*h^2+h^3).
For C!=0, reduce A modulo C times the quartics; polynomial division recovers
the unique h attaining that reduced prefix. Retain all base transports
attaining the least prefix. For C=0, reduce B modulo {A*h+h^3} instead.
Row-echelon reduction zeros the earliest pivot coordinates, giving the
lexicographically least suffix in each coset.

These steps minimize over all 6,144 F_2 or 5,832 F_3 coordinate changes.
The helper precomputes their affine suffix maps once per prefix.
Allowed changes are invertible and preserve the equation class, smoothness,
and smooth-model minimality. Equivalent equations have identical coordinate
orbits and hence equal I. Conversely, equal I supplies invertible changes
from both equations to the same normalized model. Thus I separates the
coordinate orbits counted here, including singular and nonminimal equations.
It does not classify abstract surfaces under arbitrary rational changes.

## Logical multiset hash

Parse records as (I,smooth,old_minimal), with binary flags.
The disjoint buckets are `smooth_minimal`, `smooth_nonminimal`, `delta_zero`;
`all` contains their union. Sort by the full tuple, preserving multiplicity.
For each prefix and bucket, hash this ASCII stream with LF endings:

```text
census-orbit-v1|q|canonical_prefix_integer|bucket|record_count
coefficient_digits_in_fixed_width|smooth_0_or_1|old_minimal_0_or_1
...
```

The coefficient width is exactly 32 or 27, including leading zeroes.
An empty bucket hashes its header with count zero.
Combine prefix summaries in increasing numerical prefix order:

```text
census-multiset-v1|q|bucket
canonical_prefix_integer|record_count|prefix_bucket_digest
...
```

Counts are compared separately. Production comparisons also reject repeated
I values within or between buckets and repeated canonical prefix labels.
Removing or duplicating a row changes the count; a balanced substitution
changes the sorted canonical tuples. SHA-256 collisions remain the only
hash collision limit after exact count and uniqueness checks. There is no
claimed probabilistic bound beyond the usual SHA-256 collision assumption.
Serialization and row order do not enter the mathematical comparison.

## Running the comparisons

Download the LFS fixture and run from the repository root:

```sh
make check/archive_fixtures.m check/archive_invariants.m
```

Each file has its own CI matrix job with a 30-minute budget.
The tests use `FibrationPrefixes`, `NumberOfFibrations`,
`FibrationRepresentativeKeys`, and `aInvariantsFromKey`.
The independent invariant helper imports no package-private code.
Direct polynomial substitutions check its minimum on two models per field,
including an explicitly changed representative.
The complete zero-prefix records are also reversed and reserialized as JSON.
Their hash must survive that change and detect deletion, duplication,
balanced replacement, a changed invariant, and changed classification flags.
