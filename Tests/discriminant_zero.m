// Pins the zero equation's all-place certificate and a regular cusp at infinity.
procedure test_singular_equations()
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        W := [ R | 0,0,0,0,0 ];
        kind,minimality := FibrationProperties(W);
        assert kind eq "nonregular" and minimality eq "nonminimal somewhere";
        ok,certificate := IsEverywhereMinimal(W);
        assert not ok and certificate eq <"all">;
        kind,minimality := FibrationProperties([ R | 0,0,0,0,s*t^11 ]);
        assert kind eq "regular cusp" and minimality eq "nonminimal somewhere";
    end for;
end procedure;

test_singular_equations();
