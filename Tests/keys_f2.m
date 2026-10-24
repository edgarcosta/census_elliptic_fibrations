R<s,t> := PolynomialRing(GF(2),2);
keys := [0,1,2^13-1,2^13,2^22-1,2^22,2^29-1,2^29,2^32-1];
for key in keys do
    W := FibrationFromKeyF2(R,key);
    assert W[2] eq 0;
    assert EF2BitKey(&cat[ EF2Coefficients(W[i],[2,0,6,8,12][i])
                           : i in [1,3,4,5] ]) eq key;
end for;
assert FibrationFromKeyF2(R,2^31+2^28+2^21+2^12)
       eq [R | s^2,0,s^6,s^8,s^12];
assert FibrationFromKeyF2(R,2^29+2^22+2^13+1)
       eq [R | t^2,0,t^6,t^8,t^12];
W, key := CanonicalFibrationF2([R | s*t,0,s^6+t^6,s^3*t^5,s^11*t+t^12]);
assert FibrationFromKeyF2(R,key) eq W;
W2, key2 := CanonicalFibrationF2(FibrationFromKeyF2(R,key));
assert W2 eq W and key2 eq key;
