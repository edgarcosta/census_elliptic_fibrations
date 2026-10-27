# Elliptic K3 census

This Magma package counts Weierstrass model orbits over GF(2) and GF(3),
including elliptic and quasi-elliptic K3 fibrations with a section.
It also counts equations with nonregular generic fibre and nonminimal
equations. Generic-fibre kind and minimality are reported separately.

Models are identified under changes of Weierstrass coordinates preserving
the zero section and the coefficient line bundle O(2), with automorphisms
of the base allowed. Local minimalization does not identify equations in
this census. See [the mathematics](docs/mathematics.md) for the equivalence
relation and the criteria for K3 fibrations.

## Installation

Clone or download this repository and attach its `spec` file. From the
repository root, run:

```magma
AttachSpec("spec");
```

From another directory, give `AttachSpec` the absolute path to `spec`.
No other Magma package is required. The examples below were checked with
Magma V2.30-1.

For [CHIMP](https://github.com/edgarcosta/CHIMP), add this repository as a
submodule named `census_elliptic_fibrations`. Inside the outer braces of
`CHIMP.spec`, include:

```text
census_elliptic_fibrations
{
  +spec
}
```

Attach `CHIMP.spec` as usual. This follows CHIMP's
[submodule spec entries](https://github.com/edgarcosta/CHIMP/blob/master/CHIMP.spec).

## Counts

The input ring is a binary polynomial ring over GF(2) or GF(3).
`NumberOfFibrations` returns the total, an associative array of cells,
and pairs `<prefix index, cell counts>` for the selected prefixes.

```magma
AttachSpec("spec");
R<s,t> := PolynomialRing(GF(3),2);
n,cells,rows := NumberOfFibrations(R);
print n;
print cells[<"smooth","everywhere minimal">];
print cells[<"regular cusp","everywhere minimal">];
```

```text
1307814237
1307723786
294
```

The default includes every orbit. For the elliptic K3 census, select
`GenericFibreKinds := ["smooth"]` and
`Minimalities := ["everywhere minimal"]`. For the quasi-elliptic K3
census, use `"regular cusp"` with the same minimality selection.
The [interface](docs/interface.md) describes all selections and return values.

## Representatives

Representatives are listed over one prefix at a time. A key encodes a
normalized equation; it need not be the least key in its orbit.
This example lists the quasi-elliptic K3 representatives over GF(3).

```magma
AttachSpec("spec");
R<s,t> := PolynomialRing(GF(3),2);
prefixes := FibrationPrefixes(R);
i := Index(prefixes,<R!0,R!0>);
keys := FibrationRepresentativeKeys(R,i :
    GenericFibreKinds := ["regular cusp"],
    Minimalities := ["everywhere minimal"]);
print #keys;
W := aInvariantsFromKey(R,keys[1]);
kind,minimality := FibrationProperties(W);
print kind;
print minimality;
```

```text
294
regular cusp
everywhere minimal
```

`WriteFibrationRepresentativeKeys` writes selected prefixes to binary files
with a summary and completion metadata. It holds one prefix's representatives
in memory. By default, it writes all prefixes and all six cells: smooth,
regular cusp and nonregular, each with both minimalities. Over GF(3),
that is 1,307,814,237 keys, occupying 10,462,513,896 bytes at eight bytes
per key, excluding metadata.
Choose prefixes before listing or writing representatives.
See [keys and files](docs/keys.md) for the format and a bounded writing example.

## Results

The six cells are disjoint. Only the smooth, everywhere minimal and regular
cusp, everywhere minimal cells describe K3 fibrations with a section.

| Generic-fibre kind | Minimality | GF(2) | GF(3) |
| --- | --- | ---: | ---: |
| smooth | everywhere minimal | 1,397,598 | 1,307,723,786 |
| smooth | nonminimal somewhere | 4,134 | 90,113 |
| regular cusp | everywhere minimal | 680 | 294 |
| regular cusp | nonminimal somewhere | 38 | 20 |
| nonregular | everywhere minimal | 7 | 12 |
| nonregular | nonminimal somewhere | 11 | 12 |
| All | | 1,402,468 | 1,307,814,237 |

These counts are reproduced by `NumberOfFibrations` over both fields.
The [mathematics](docs/mathematics.md#counts-and-independent-checks) describes
the independent checks and their scope.

## Tests and examples

Thanks for helping check the census. From the repository root, run either
fast-suite command:

```sh
magma -b exitsignal:="" run_tests.m
make check
```

Both run `Tests/*.m` and `examples/*.m` and return a nonzero status on failure.
The examples include counts over both fields and a bounded representative
calculation. The CI matrix runs these same files.
The tests need Magma, a POSIX shell, standard shell utilities, and `python3`
for the reference launcher.

To run one file with the runner's `filename` parameter:

```sh
magma -b filename:=Tests/keys_f3.m exitsignal:="" run_tests.m
```

Use `target` to select paths containing a string. The `debug` parameter opens
the debugger on an error. The longer comparisons have separate Makefile targets,
described in [the reference guide](reference/README.md).

## Documentation

- [Interface](docs/interface.md): intrinsics, selections, and certificates.
- [Mathematics](docs/mathematics.md): normal forms, group actions, counting, and classification.
- [Keys and files](docs/keys.md): normalization, coefficient encodings and saved representatives.
- [Reference calculations](reference/README.md): independent scan, Burnside, and local checks.

Author: Edgar Costa
