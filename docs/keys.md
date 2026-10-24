# Keys and representative files

A key records the coefficients of a normalized equation over GF(2) or GF(3).
Keep its characteristic with it. The same integer has different meanings
over the two fields. A key contains no prefix index, property or orbit
identifier, and is not canonical under changes of coordinates.

The [interface](interface.md) gives the conversion and writing intrinsics;
the [mathematics](mathematics.md) defines the model equivalence relation.

## Normalization

For W = [a1,a2,a3,a4,a6], translation x_old = x+a2 over GF(2) gives

```text
N(W) = [a1,0,a3+a1*a2,a4+a2^2,a6+a2*a4].
```

Over GF(3), completing the square gives

```text
C = a2+a1^2/4,
A = a4+a1*a3/2,
B = a6+a3^2/4,
N(W) = [0,C,0,A,B].
```

Decoding FibrationKey(W) gives N(W). Normalization makes no further
choice among the models in an orbit.

## Integer keys

Key encoding version 1 uses the coefficient basis

```text
s^d, s^(d-1)*t, ..., s*t^(d-1), t^d
```

for a degree-d form, including all zero coefficients. Here s = R.1 and
t = R.2. Concatenate the coefficient blocks in the following order and
read the resulting digits in base p, most significant digit first.

| Field | Blocks | Digits | Range |
| --- | --- | ---: | --- |
| GF(2) | a1, a3, a4, a6 | 3+7+9+13 = 32 | 0 <= key < 2^32 |
| GF(3) | C, A, B | 5+9+13 = 27 | 0 <= key < 3^27 |

Every integer in the range encodes a normalized equation. The printed
variable names do not affect the encoding; exchanging the generators does.
These examples pin the first coefficient of each block:

```magma
AttachSpec("spec");
R<s,t> := PolynomialRing(GF(2),2);
key := 2^31+2^28+2^21+2^12;
W := aInvariantsFromKey(R,key);
print W;
print FibrationKey(W) eq key;
R<s,t> := PolynomialRing(GF(3),2);
key := 3^26+3^21+3^12;
W := aInvariantsFromKey(R,key);
print W;
print FibrationKey(W) eq key;
```

```text
[ s^2, 0, s^6, s^8, s^12 ]
true
[ 0, s^4, 0, s^8, s^12 ]
true
```

## Representative files

For a stem and prefix index i, the writer creates:

| Path | Contents |
| --- | --- |
| stem.csv | One row per selected cell per prefix |
| stem_iiii.bin | The selected representative keys, eight bytes per key |
| stem_iiii.stats | Prefix forms, selections, counts, and completion marker |

The filename index is padded to four digits. The containing directory
must exist. All output paths are checked before writing: existing binary,
sidecar, or summary files cause an error. Concurrent calls need distinct stems.
There is no resume operation.

This bounded example writes only the nonregular equation in the zero prefix
over GF(3). Create the ignored `artifacts/` directory first, then run from
the repository root with a fresh stem:

```magma
AttachSpec("spec");
R<s,t> := PolynomialRing(GF(3),2);
i := Index(FibrationPrefixes(R),<R!0,R!0>);
n := WriteFibrationRepresentativeKeys(R,"artifacts/documented-zero" :
    PrefixIndices := [i], GenericFibreKinds := ["nonregular"]);
print n;
prefix_stem := Sprintf("artifacts/documented-zero_%o",IntegerToString(10000+i)[2..5]);
F := Open(prefix_stem cat ".bin","rb");
print ReadBytes(F);
delete F;
print Read(prefix_stem cat ".stats");
```

In the current prefix table i = 1. The output is

```text
1
[ 0, 0, 0, 0, 0, 0, 0, 0 ]
file_format_version=1
key_encoding_version=1
characteristic=3
prefix_index=1
prefix_form_1_degree=4
prefix_form_1_key=0
prefix_form_2_degree=8
prefix_form_2_key=0
generic_fibre_kinds=nonregular
minimalities=everywhere minimal|nonminimal somewhere
number_of_keys=1
cell=nonregular|everywhere minimal|0
cell=nonregular|nonminimal somewhere|1
complete=1
```

### Binary bytes

The binary file has no header. For integer key k, byte j is
(k div 256^j) mod 256 for j = 0,...,7: unsigned, little-endian.
A length not divisible by eight is incomplete or invalid. Divisibility
alone does not establish completion. The sidecar supplies the characteristic
and key-encoding version. A prefix file combines all selected cells;
decode a key and call FibrationProperties to recover its properties.

### Text formats

The current file format is version 1. Text is ASCII with LF line endings
and a final LF. Integers are decimal without leading zeros, except for
the filename index. Fields contain no quoting or escaping.

The sidecar lines have the order shown in the example. There is one
`cell=kind|minimality|count` line for every selected cell, including zeros.
Selections and cell lines follow the property order in the
[interface](interface.md#selections), regardless of the order of the options.
Empty option sequences expand to the full value lists.
The form degrees and base-p keys specify the prefix without variable names.
`number_of_keys` equals the sum of the cell counts.

The CSV header is

```text
file_format_version,key_encoding_version,characteristic,prefix_index,prefix_form_1_degree,prefix_form_1_key,prefix_form_2_degree,prefix_form_2_key,generic_fibre_kind,minimality,number_of_keys
```

Rows follow the caller's prefix order, then generic-fibre kind and minimality
in property order. The prefix metadata agrees with the sidecar.
The last field counts that row's cell. Even a prefix with no matching
orbits gets an empty binary file, its selected zero-count rows, and a
completed sidecar.

### Completion

For each prefix the writer flushes and closes the binary file, appends
and flushes its CSV rows, then writes and closes the sidecar with
`complete=1` last. A completed prefix has a fully parsed sidecar ending
in that marker, matching CSV rows, and binary length 8*number_of_keys.
An interrupted prefix has no completion guarantee. Ignore summary rows
without a completed sidecar. An output error raises an error; it does
not return a successful partial count.

Representative order and the representative chosen within an orbit may
change between package versions. Keep the prefix forms and both format
versions when transferring files; a bare prefix index is not permanent.
The plain reference routines have their own historical output format.
