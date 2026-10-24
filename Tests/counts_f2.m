R<s,t> := PolynomialRing(GF(2),2);
prefixes := FibrationOrbitPrefixesF2G(R);
assert #prefixes eq 88;
counts, totals := EF_CensusCounts(prefixes, FibrationFibreF2G);
assert #counts eq #prefixes;
assert totals eq <1402468,1401732,1397598>;
