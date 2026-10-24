R<s,t> := PolynomialRing(GF(3),2);
keys := [0,1,3^13-1,3^13,3^22-1,3^22,3^27-1];
for key in keys do
    assert EF3_ModelKey(FibrationFromKeyF3(R,key)) eq key;
end for;
assert FibrationFromKeyF3(R,3^26+3^21+3^12) eq [R | s^4,s^8,s^12];
assert FibrationFromKeyF3(R,3^22+3^13+1) eq [R | t^4,t^8,t^12];
W := [R | s^4+2*s*t^3,2*s^7*t+t^8,s^12+2*s^5*t^7+t^12];
assert FibrationFromKeyF3(R,EF3_ModelKey(W)) eq W;
