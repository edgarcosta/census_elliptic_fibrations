import "tests/support/_support.m": expect_error;

// Pins coefficient order, block boundaries and normalization before encoding.
procedure test_key_encoding_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    assert aInvariantsFromKey(R, 3^26+3^21+3^12) eq [R | 0, s^4, 0, s^8, s^12];
    assert aInvariantsFromKey(R, 3^22+3^13+1) eq [R | 0, t^4, 0, t^8, t^12];
    P<z> := PolynomialRing(GF(3));
    degrees := [2, 4, 6, 8, 12];
    indices := [2, 4, 5];
    for i in indices do
        for j in [0..degrees[i]] do
            W := [R | 0, 0, 0, 0, 0];
            W[i] := (3-1)*s^(degrees[i]-j)*t^j;
            key := (3-1)*3^(&+([0] cat [degrees[k]+1 : k in indices | k gt i])+degrees[i]-j);
            assert aInvariantsFromKey(R, key) eq W;
            assert Evaluate(W[i], [z, P!1]) eq (3-1)*z^(degrees[i]-j);
            assert Evaluate(W[i], [P!1, z]) eq (3-1)*z^j;
        end for;
    end for;
    keys := [0, 1] cat &cat[[3^d-1, 3^d] : d in [13, 22, 27]];
    Prune(~keys);
    assert #keys gt 0;
    for key in keys do
        assert FibrationKey(aInvariantsFromKey(R, key)) eq key;
    end for;
    W := [R | s^2, s^3*t, s^5*t, s^7*t, s^11*t];
    assert aInvariantsFromKey(R, FibrationKey(W)) eq [R | 0, W[2]+W[1]^2/4, 0, W[4]+W[1]*W[3]/2, W[5]+W[3]^2/4];
    expect_error(procedure() x := aInvariantsFromKey(R, -1); end procedure, "key is outside");
    expect_error(procedure() x := aInvariantsFromKey(R, 3^27); end procedure, "key is outside");
    expect_error(procedure() x := FibrationKey([R | s, 0, 0, 0, 0]); end procedure, "five homogeneous coefficients");
    expect_error(procedure() x := FibrationKey([R | s^2+1, 0, 0, 0, 0]); end procedure, "five homogeneous coefficients");
    expect_error(procedure() x := FibrationKey([R | 0, 0, 0]); end procedure, "five homogeneous coefficients");
    expect_error(procedure() x := FibrationKey([0, 0, 0, 0, 0]); end procedure, "five homogeneous coefficients");
    S := PolynomialRing(GF(3^2), 2);
    expect_error(procedure() x := aInvariantsFromKey(S, 0); end procedure, "binary polynomial ring");
    S := PolynomialRing(GF(3), 3);
    expect_error(procedure() x := aInvariantsFromKey(S, 0); end procedure, "binary polynomial ring");
end procedure;

test_key_encoding_f3();

// Pins both charts and weighted coefficients without normalization or minimization.
// Built-in forms share these inputs to check overload registration and dispatch.
procedure test_model_conversions_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    F<u> := FunctionField(GF(3));
    // Rank-one function fields can have multivariate polynomial numerators.
    G := FunctionField(GF(3), 1);
    fixtures := [
        <[F | 0, 1, 0, 2, 0],
         [R | 0, t^4, 0, 2*t^8, 0],
         [F | 0, u^4, 0, 2*u^8, 0]>,
        <[F | 2*u^2+u+1, u^4+2*u, 2*u^6+u^2+1, u^8+2*u^3,
                2*u^12+u+1],
         [R | 2*s^2+s*t+t^2, s^4+2*s*t^3, 2*s^6+s^2*t^4+t^6,
                s^8+2*s^3*t^5, 2*s^12+s*t^11+t^12],
         [F | 2+u+u^2, 1+2*u^3, 2+u^4+u^6, 1+2*u^5,
                2+u^11+u^12]>,
        <[F | 0, 0, 0, u^8, u],
         [R | 0, 0, 0, s^8, s*t^11],
         [F | 0, 0, 0, 1, u^11]>];
    assert #fixtures gt 0;
    for fixture in fixtures do
        A, W, reciprocal := Explode(fixture);
        E := EllipticCurve(A);
        assert aInvariants(E) eq A;
        assert aInvariants(E, R) eq W;
        assert aInvariants(EllipticCurve(W, F)) eq A;
        assert aInvariants(EllipticCurve(W, G), R) eq W;
        Einf := EllipticCurve(W, F : Infinity := true);
        assert aInvariants(Einf) eq reciprocal;
        assert aInvariants(EllipticCurve(reciprocal)) eq reciprocal;
    end for;
    assert aInvariants(EllipticCurve(fixtures[1][3]), R)
        eq [R | 0, s^4, 0, 2*s^8, 0];
end procedure;

test_model_conversions_f3();

procedure test_conversion_degree_bounds_f3()
    R := PolynomialRing(GF(3), 2);
    F<u> := FunctionField(GF(3));
    degrees := [2, 4, 6, 8, 12];
    for i in [1..5] do
        A := [F | 0, 1, 0, 2, 0];
        A[i] := u^(degrees[i]+1);
        E := EllipticCurve(A);
        expect_error(procedure() x := aInvariants(E, R); end procedure,
            "coefficients must be polynomials of degrees at most 2,4,6,8,12");
    end for;
end procedure;

test_conversion_degree_bounds_f3();

procedure test_conversion_poles_f3()
    R := PolynomialRing(GF(3), 2);
    F<u> := FunctionField(GF(3));
    for i in [1..5] do
        A := [F | 0, 1, 0, 2, 0];
        A[i] := 1/u;
        E := EllipticCurve(A);
        expect_error(procedure() x := aInvariants(E, R); end procedure,
            "coefficients must be polynomials of degrees at most 2,4,6,8,12");
    end for;
end procedure;

test_conversion_poles_f3();

procedure test_conversion_model_shape_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    F := FunctionField(GF(3));
    for a6 in [s, s^12+s] do
        W := [R | 0, t^4, 0, 2*t^8, 0];
        W[5] := a6;
        expect_error(procedure() E := EllipticCurve(W, F); end procedure,
            "five homogeneous coefficients");
    end for;
    expect_error(procedure() E := EllipticCurve([R | 0, 0, 0], F);
        end procedure, "five homogeneous coefficients");
    expect_error(procedure() E := EllipticCurve([0, 0, 0, 0, 0], F);
        end procedure, "five homogeneous coefficients");
    S := PolynomialRing(GF(3), 3);
    expect_error(procedure() E := EllipticCurve([S | 0, 0, 0, 0, 0], F);
        end procedure, "five homogeneous coefficients");
    E := EllipticCurve([F | 0, 1, 0, 2, 0]);
    expect_error(procedure() W := aInvariants(E, S); end procedure,
        "binary polynomial ring");
end procedure;

test_conversion_model_shape_f3();

procedure test_conversion_constant_fields_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    H := FunctionField(GF(3^2));
    S<v, w> := PolynomialRing(GF(3^2), 2);
    W := [R | 0, t^4, 0, 2*t^8, 0];
    E := EllipticCurve([H | 0, 1, 0, 2, 0]);
    expect_error(procedure() C := EllipticCurve(W, H); end procedure,
        "same constant field");
    expect_error(procedure() A := aInvariants(E, R); end procedure,
        "same constant field");
    V := [S | Evaluate(a, [v, w]) : a in W];
    expect_error(procedure() C := EllipticCurve(V, H); end procedure,
        "five homogeneous coefficients");
    expect_error(procedure() A := aInvariants(E, S); end procedure,
        "binary polynomial ring");
end procedure;

test_conversion_constant_fields_f3();

procedure test_conversion_singular_models_f3()
    R<s, t> := PolynomialRing(GF(3), 2);
    F := FunctionField(GF(3));
    models := [[R | 0, 0, 0, 0, s*t^11], [R | 0, 0, 0, 0, 0],
        [R | 0, t^4, 0, 0, 0],
        // Completing the square gives y^2 = x*(x-1)^2 in the finite chart.
        [R | t^2, 0, 0, t^8, 0]];
    assert #models gt 0;
    for W in models do
        expect_error(procedure() E := EllipticCurve(W, F); end procedure,
            "W must have nonzero discriminant");
    end for;
end procedure;

test_conversion_singular_models_f3();
