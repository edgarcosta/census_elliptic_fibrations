import "Tests/_support.m": expect_error;

// Pins coefficient order, block boundaries and normalization before encoding.
procedure test_key_encoding_f3()
    R<s,t> := PolynomialRing(GF(3),2);
    assert aInvariantsFromKey(R,3^26+3^21+3^12) eq [ R | 0,s^4,0,s^8,s^12 ];
    assert aInvariantsFromKey(R,3^22+3^13+1) eq [ R | 0,t^4,0,t^8,t^12 ];
    keys := [0,1] cat &cat[ [3^d-1,3^d] : d in [13, 22, 27] ];
    Prune(~keys);
    assert #keys gt 0;
    for key in keys do
        assert FibrationKey(aInvariantsFromKey(R,key)) eq key;
    end for;
    W := [ R | s^2,s^3*t,s^5*t,s^7*t,s^11*t ];
    assert aInvariantsFromKey(R,FibrationKey(W)) eq [ R | 0,W[2]+W[1]^2/4,0,W[4]+W[1]*W[3]/2,W[5]+W[3]^2/4 ];
    expect_error(procedure() x := aInvariantsFromKey(R,-1); end procedure,"key is outside");
    expect_error(procedure() x := aInvariantsFromKey(R,3^27); end procedure,"key is outside");
    expect_error(procedure() x := FibrationKey([ R | s,0,0,0,0 ]); end procedure,"five homogeneous coefficients");
    expect_error(procedure() x := FibrationKey([ R | 0,0,0 ]); end procedure,"five homogeneous coefficients");
    expect_error(procedure() x := FibrationKey([ 0,0,0,0,0 ]); end procedure,"five homogeneous coefficients");
    S := PolynomialRing(GF(3^2),2);
    expect_error(procedure() x := aInvariantsFromKey(S,0); end procedure,"binary polynomial ring");
    S := PolynomialRing(GF(3),3);
    expect_error(procedure() x := aInvariantsFromKey(S,0); end procedure,"binary polynomial ring");
end procedure;

test_key_encoding_f3();
