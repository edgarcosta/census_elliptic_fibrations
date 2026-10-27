// Pins finite, empty and all-place certificates, including infinity.
procedure test_minimality_certificates_f2()
    R<s,t> := PolynomialRing(GF(2),2);
    ok,certificate := IsEverywhereMinimal([ R | 0,0,0,0,s*t^11 ]);
    assert not ok and certificate eq <"finite",[ R | t ]>;
    ok,certificate := IsEverywhereMinimal([ R | 0,0,0,s^6*t^2,0 ]);
    assert not ok and certificate eq <"finite",[ R | s ]>;
    ok,certificate := IsEverywhereMinimal([ R | s*t,s^4+t^4,0,0,0 ]);
    assert ok and certificate eq <"finite",[ R | ]>;
    ok,certificate := IsEverywhereMinimal([ R | s*t,0,0,0,s^6*t^6 ]);
    assert not ok and certificate eq <"finite",[ R | t,s ]>;
    pi := s^2+s*t+t^2;
    ok,certificate := IsEverywhereMinimal([ R | pi,0,0,0,pi^6 ]);
    assert not ok and certificate eq <"finite",[ R | pi ]>;
    for W in [[ R | 0,0,0,0,0 ],[ R | 0,0,0,s^8,t^12 ]] do
        ok,certificate := IsEverywhereMinimal(W);
        assert not ok and certificate eq <"all">;
    end for;
end procedure;

test_minimality_certificates_f2();
