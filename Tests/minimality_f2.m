R<s,t> := PolynomialRing(GF(2),2);
W := [R | s*t,0,0,0,s^6*(s^6+t^6)];
assert FibrationDiscriminantF2(W) ne 0;
ok, bad := IsEverywhereMinimalF2(W);
assert not ok and bad eq [R | s];
ok, bad := IsEverywhereMinimalF2([R | Evaluate(a,[t,s]) : a in W]);
assert not ok and bad eq [R | t];
ok, bad := IsEverywhereMinimalF2([R | s*t,0,0,0,s^12+t^12]);
assert ok and IsEmpty(bad);
