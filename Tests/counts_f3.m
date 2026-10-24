R<s,t> := PolynomialRing(GF(3),2);
prefixes := FibrationOrbitPrefixesF3G(R);
assert #prefixes eq 1860;
counts, totals := EF_CensusCounts(prefixes, FibrationFibreF3G);
assert #counts eq #prefixes;
assert totals eq <1307814237,1307813899,1307723786>;
