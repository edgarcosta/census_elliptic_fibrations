# Interface

The public intrinsics count model orbits, list representatives, and test
individual equations. Attach `spec` before calling them.
The [README](../README.md) gives installation and initial examples.

## Rings and models

`R` is a polynomial ring in two variables over exactly GF(2) or GF(3).
Extension fields and rings of other ranks are rejected.
Write `s = R.1` and `t = R.2`; their order determines the coefficient encoding.

A model is a sequence `W = [a1,a2,a3,a4,a6]` in `R`, with homogeneous
degrees `[2,4,6,8,12]`, for the equation

```text
y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6.
```

Zero coefficients are allowed. Use a typed sequence such as
`[R | 0,0,0,0,0]` for the zero equation.
An integer sequence or a characteristic-3 triple `[C,A,B]` is not a public model.
Invalid arguments raise a requirement error with a message for the caller.

## Selections

Counts, listings and file output accept these case-sensitive property values:

| Option | Values, in table order |
| --- | --- |
| `GenericFibreKinds` | `"smooth"`, `"regular cusp"`, `"nonregular"` |
| `Minimalities` | `"everywhere minimal"`, `"nonminimal somewhere"` |

Each option defaults to `[]`, which selects every value.
A nonempty sequence selects its distinct listed values; repeats and unknown
values are rejected. The selections are intersected.
A regular cusp is regular over the function field with geometrically
cuspidal generic fibre.

`PrefixIndices := []` selects every prefix for counts and writing.
A nonempty sequence gives distinct valid indices in the order to process them.
Listings take one prefix index as a positional argument.
All selections are checked before counting or creating output files.

## Prefixes and counts

> **FibrationPrefixes** (R :: *RngMPol*) -> *SeqEnum*

The ordered pairs `<a1,a3>` over GF(2), or `<C,A>` after completing the
square over GF(3). The tables have 88 and 1,860 entries, respectively.
A prefix is an intermediate orbit of coefficient pairs, usually containing
many full model orbits. Indices belong to the installed package's table;
they are not permanent identifiers and need not match the reference indices.

> **NumberOfFibrations** (R :: *RngMPol* : PrefixIndices := [],
> GenericFibreKinds := [], Minimalities := []) -> *RngIntElt*, *Assoc*, *SeqEnum*

The selected orbit count `n`, the cell counts `cells`, and the prefix rows
`rows`. A cell has key `<kind,minimality>`. Only selected cells are present,
including cells with count zero. Each row is `<i,prefix_cells>` with the
same selected cell keys. The total is the sum of all cells, and each cell
is the sum of that cell over the rows.

This example selects quasi-elliptic K3 fibrations over both fields:

```magma
AttachSpec("spec");
for p in [2,3] do
    R := PolynomialRing(GF(p),2);
    n,cells,rows := NumberOfFibrations(R :
        GenericFibreKinds := ["regular cusp"],
        Minimalities := ["everywhere minimal"]);
    printf "GF(%o): %o prefixes, %o quasi-elliptic K3 orbits\n",
        p,#FibrationPrefixes(R),n;
end for;
```

```text
GF(2): 88 prefixes, 680 quasi-elliptic K3 orbits
GF(3): 1860 prefixes, 294 quasi-elliptic K3 orbits
```

## Representatives and keys

> **FibrationRepresentativeKeys** (R :: *RngMPol*, i :: *RngIntElt* :
> GenericFibreKinds := [], Minimalities := []) -> *SeqEnum*

One integer key per selected orbit over prefix `i`.
An empty result is a typed integer sequence.
The length agrees with the selected count for that prefix.
The result is held in memory; representative choice and order can vary
between package versions. See the [README example](../README.md#representatives).

> **aInvariantsFromKey** (R :: *RngMPol*, key :: *RngIntElt*) -> *SeqEnum*

The normalized five-coefficient model encoded by `key`.
Every integer in the [encoding range](keys.md#integer-keys) is accepted;
decoding makes no assertion about generic-fibre kind or minimality.

> **FibrationKey** (W :: *SeqEnum*) -> *RngIntElt*

The coefficient key after normalization. It records the normalized equation,
without searching its orbit. Decoding recovers that normalization of `W`.
Encoding a decoded valid key recovers the same integer.

> **WriteFibrationRepresentativeKeys** (R :: *RngMPol*, stem :: *MonStgElt* :
> PrefixIndices := [], GenericFibreKinds := [], Minimalities := []) -> *RngIntElt*

The number of keys written to new files with the given stem.
The containing directory must exist. Existing output paths are rejected
before any file is created. The writer processes one prefix at a time;
it keeps that prefix's representatives in memory.
See [keys and files](keys.md#representative-files) for the completion rules.

## Properties and minimality

> **FibrationProperties** (W :: *SeqEnum*) -> *MonStgElt*, *MonStgElt*

The generic-fibre kind followed by the minimality string, using the values
above. These are separate properties. A nonregular equation may be
everywhere minimal. A regular cusp may be nonminimal somewhere.

> **IsEverywhereMinimal** (W :: *SeqEnum*) -> *BoolElt*, *Tup*

Whether the equation is minimal at every closed place of the base,
followed by a certificate. Nonregular equations are accepted.

| Boolean | Certificate | Meaning |
| --- | --- | --- |
| `true` | `<"finite",[R \| ]>` | Everywhere minimal |
| `false` | `<"finite",places>` | Nonminimal exactly at the nonempty list of places |
| `false` | `<"all">` | Nonminimal at every closed place |

Finite places are represented by homogeneous irreducible forms, with first
nonzero coefficient equal to one. They are ordered by degree and then by
their coefficient key. The place at infinity in the chart `a_i(s,1)` is `t`.
The certificate always includes infinity when it is bad.

Here are a regular cusp nonminimal only at infinity and the zero equation:

```magma
AttachSpec("spec");
R<s,t> := PolynomialRing(GF(3),2);
W := [R | 0,0,0,0,s*t^11];
kind,minimality := FibrationProperties(W);
print kind;
print minimality;
minimal,certificate := IsEverywhereMinimal(W);
print minimal;
print certificate;
minimal,certificate := IsEverywhereMinimal([R | 0,0,0,0,0]);
print minimal;
print certificate;
```

```text
regular cusp
nonminimal somewhere
false
<"finite", [
    t
]>
false
<"all">
```

The [mathematics](mathematics.md#local-minimality) gives the divisibility
criterion and explains why finitely many candidate places suffice.

## Progress messages

`SetVerbose("FibrationCensus",1)` reports completed or written prefixes.
Level 2 also reports fibre dimension, group order and CPU time for counts.
The default is silent.
