// The zero prefix contains every quasi-elliptic K3 orbit in each characteristic.
procedure example_representatives()
    for p in [2,3] do
        R := PolynomialRing(GF(p),2);
        i := Index(FibrationPrefixes(R),<R!0,R!0>);
        keys := FibrationRepresentativeKeys(R,i :
            GenericFibreKinds := ["regular cusp"],
            Minimalities := ["everywhere minimal"]);
        assert #keys eq (p eq 2 select 680 else 294);
        W := aInvariantsFromKey(R,keys[1]);
        kind,minimality := FibrationProperties(W);
        assert kind eq "regular cusp" and minimality eq "everywhere minimal";
        assert FibrationKey(W) eq keys[1];
    end for;
end procedure;

example_representatives();
