R<s,t> := PolynomialRing(GF(3),2);
W := [R | s^2*t^2,0,s^6*(s^6+t^6)];
assert FibrationDiscriminantF3(W) ne 0;
ok, bad := IsEverywhereMinimalF3(W);
assert not ok and bad eq [R | s];
ok, bad := IsEverywhereMinimalF3([R | Evaluate(a,[t,s]) : a in W]);
assert not ok and bad eq [R | t];

// The affine model y^2 = x^3 + u^4*x + u is minimal at u=0 despite v(Delta)=12.
// Its degree-(4,8,12) homogenization is non-minimal only at infinity.
W := [R | 0,s^4*t^4,s*t^11];
assert FibrationDiscriminantF3(W) eq -s^12*t^12;
ok, bad := IsEverywhereMinimalF3(W);
assert not ok and bad eq [R | t];
ok, bad := IsEverywhereMinimalF3([R | 0,s^4*t^4,s*t^11+s^11*t]);
assert ok and IsEmpty(bad);
