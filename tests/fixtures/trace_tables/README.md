# Prototype trace tables

These fixtures let you compare the exact Magma computation with the existing
ordinary tables. The F_3 fixture values were built by the prototype's complex
Fourier transform. They are comparison data; Magma's `TraceOfFrobenius` provides the
independent check through degree 4.

The source is `artifacts/mw/lfun-counter/field.cpp` with `field.hpp` in the
native-stage-1 prototype snapshot. `sources.sha256` records both source hashes.
An exporter constructed `Field(p,d)` for p in [2,3] and d in [1..12], then wrote
`ordinary[0..p^d-2]`. It omitted the two unused trailing slots.
The polynomial coefficients were copied from that source's `polys2` and `polys3`.

Each ASCII file has three lines, each ending with LF:

1. The characteristic and extension degree, separated by one space.
2. The monic field polynomial's coefficients, constant term first.
3. The signed integer traces at alpha^e, for e from 0 through p^d-2.

Here alpha is the polynomial residue class, including alpha=2 in degree 1 over
F_3. The curves are y^2+x*y=x^3+a over F_2 and y^2=x^3+x^2+a over F_3.
Entries are separated by one space, with no leading or trailing spaces.
Parameter zero is omitted; its nodal trace is 1 in both characteristics.

`tables.sha256` hashes the complete files, including the field definitions.
The test serializes the computed table in this format and checks each hash.
All table files use Git LFS.
