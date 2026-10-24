// Verifies a bounded quasi-elliptic listing without collecting a whole census.
procedure example_representatives()
    R := PolynomialRing(GF(3),2);
    keys := FibrationRepresentativeKeys(R,1 : GenericFibreKinds := ["regular cusp"],
        Minimalities := ["everywhere minimal"]);
    assert #keys eq 294;
    kind,minimality := FibrationProperties(aInvariantsFromKey(R,keys[1]));
    assert kind eq "regular cusp" and minimality eq "everywhere minimal";
end procedure;

example_representatives();
