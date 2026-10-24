for p in [2,3] do
    R<s,t> := PolynomialRing(GF(p),2);
    prefixes := p eq 2 select FibrationOrbitPrefixesF2G(R)
                      else FibrationOrbitPrefixesF3G(R);
    P := prefixes[2];
    F := p eq 2 select FibrationFibreF2G(P) else FibrationFibreF3G(P);
    keys := EF_FibreRepresentatives(F);
    assert #keys eq EF_FibreCounts(F)[3];
    assert #SequenceToSet(keys) eq #keys;
    for key in keys do
        assert F`Key(F`FromKey(key)) eq key;
        model := p eq 2 select FibrationFromKeyF2(R,key) else FibrationFromKeyF3(R,key);
        assert p eq 2 select IsEverywhereMinimalF2(model) else IsEverywhereMinimalF3(model);
    end for;
end for;

k := GF(3);
d := 11;
W := VectorSpace(k,d);
b := W.d;
E := sub< GL(d+1,k) | EF_AffineMatrix(IdentityMatrix(k,d),b) >;
F := rec< EF_Fibre | Dimension := d, Group := E,
    Singular := <W!0,sub< W | [ W.i : i in [2..d] ] >>,
    NonMinimal := [ <W.2,sub< W | [ W.i : i in [1..d] | i ne 2 ] >> ],
    Key := func< v | 3^12+EF_PrimeFieldKey(v,3) > >;
keys := EF_FibreRepresentatives(F);
assert #keys eq 4*3^8;
assert SequenceToSet(keys) eq { 3^12+n : n in [0..3^d-1] |
    n div 3^(d-1) ne 0 and (n div 3^(d-2)) mod 3 ne 1 and n mod 3 eq 0 };
