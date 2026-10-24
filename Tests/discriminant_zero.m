R<s,t> := PolynomialRing(GF(2),2);
for W in [[R | 0,0,0,0,0], [R | 0,0,0,s^8,s^12]] do
    assert FibrationDiscriminantF2(W) eq 0;
    ok, bad := IsEverywhereMinimalF2(W);
    assert not ok and IsEmpty(bad);
end for;
R<s,t> := PolynomialRing(GF(3),2);
for W in [[R | 0,0,0], [R | 0,0,s^12]] do
    assert FibrationDiscriminantF3(W) eq 0;
    ok, bad := IsEverywhereMinimalF3(W);
    assert not ok and IsEmpty(bad);
end for;
