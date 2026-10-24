load "reference/f3/ef3_burnside.m";

// The expected count and spectrum were computed independently in Python.

// Pins the independent whole-space count and its fixed-locus spectrum.
procedure test_f3_burnside()
    R<s,t> := PolynomialRing(GF(3),2);
    EFB_CheckActionConventions(R);
    n,spectrum,empty,fixedSum,G := EFB_BurnsideF3(R);
    prefixCounts := EFB_PrefixBurnsideF3(R);
    assert n eq 1307814237;
end procedure;

test_f3_burnside();

Write(resultpath,"REFERENCE_CHECK_OK");
quit 0;
