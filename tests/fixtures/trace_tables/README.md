# Trace table digests over F_2 and F_3

`tables.sha256` covers all 24 tables for characteristics 2 and 3 through
extension degree 12. The test serializes each computed table as three ASCII
lines, each ending with LF:

1. The characteristic and extension degree, separated by one space.
2. The monic field polynomial's coefficients, constant term first.
3. The signed integer traces at alpha^e, for e from 0 through p^d-2.

Here alpha is the polynomial residue class, including alpha=2 in degree 1 over
F_3. The curves are y^2+x*y=x^3+a over F_2 and y^2=x^3+x^2+a over F_3.
Entries are separated by one space, with no leading or trailing spaces.
Parameter zero is omitted; its nodal trace is 1 in both characteristics.

The digests cover the field definitions and every trace.
