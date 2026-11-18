// Pins finite, empty and all-place certificates, including infinity.
procedure test_minimality_certificates_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    ok, certificate := IsEverywhereMinimal([R | 0, 0, 0, 0, s*t^11]);
    assert not ok and certificate eq <"finite", [R | t]>;
    ok, certificate := IsEverywhereMinimal([R | 0, s^4, 0, 0, 0]);
    assert not ok and certificate eq <"finite", [R | s]>;
    ok, certificate := IsEverywhereMinimal([R | 0, s*t*(s^2+t^2), 0, 0, 0]);
    assert ok and certificate eq <"finite", [R |]>;
    ok, certificate := IsEverywhereMinimal([R | 0, s^2*t^2, 0, 0, s^6*t^6]);
    assert not ok and certificate eq <"finite", [R | t, s]>;
    pi := s^2+t^2;
    ok, certificate := IsEverywhereMinimal([R | 0, pi^2, 0, 0, pi^6]);
    assert not ok and certificate eq <"finite", [R | pi]>;
    for W in [[R | 0, 0, 0, 0, 0], [R | 0, 0, 0, 0, s^12]] do
        ok, certificate := IsEverywhereMinimal(W);
        assert not ok and certificate eq <"all">;
    end for;
end procedure;

test_minimality_certificates_f3();

// Scaling x = pi^2*X, y = pi^3*Y gives a nodal model minimal away from pi.
procedure test_singular_quadratic_place_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    pi := s^2+t^2;
    W := [R | 0, pi^2, 0, 0, 0];
    ok, certificate := IsEverywhereMinimal(W);
    assert not ok and certificate eq <"finite", [R | pi]>;
    kind, minimality := FibrationProperties(W);
    assert kind eq "nonregular" and minimality eq "nonminimal somewhere";
end procedure;

test_singular_quadratic_place_f3();
