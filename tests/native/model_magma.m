SetQuitOnError(true);
AttachSpec("spec");
import "tests/_support.m": expect_error;
import "tests/_native_model_samples.m": model_samples;

// Fixed forms and base-p sums pin coefficient order independently of round trips.
procedure test_model_keys_and_cells()
    for p in [2, 3] do
        R<s, t> := PolynomialRing(GF(p), 2);
        samples := model_samples(R);
        expected := [<kind, minimality> : minimality in
            ["everywhere minimal", "nonminimal somewhere"], kind in
            ["smooth", "regular cusp", "nonregular"]];
        assert #samples eq 6 and #Seqset(expected) eq 6;
        for i in [1..#samples] do
            key, W := Explode(samples[i]);
            decoded := aInvariantsFromKey(R, key);
            assert decoded eq W and FibrationKey(decoded) eq key;
            kind, minimality := FibrationProperties(decoded);
            assert <kind, minimality> eq expected[i];
            swapped := [R | Evaluate(a, [t, s]) : a in decoded];
            kind, minimality := FibrationProperties(swapped);
            assert <kind, minimality> eq expected[i];
        end for;
        ok, bad := IsEverywhereMinimal(samples[4][2]);
        assert not ok and bad eq <"finite", [R | t]>;
        ok, bad := IsEverywhereMinimal([R | Evaluate(a, [t, s]) : a in samples[4][2]]);
        assert not ok and bad eq <"finite", [R | s]>;
    end for;
end procedure;

// The charts glue by u=1/t, x_infinity=u^4*x, y_infinity=u^6*y.
procedure test_model_charts()
    for p in [2, 3] do
        R<s, t> := PolynomialRing(GF(p), 2);
        P<z> := PolynomialRing(GF(p));
        K<u> := FunctionField(GF(p));
        S<x, y> := PolynomialRing(K, 2);
        degrees := [2, 4, 6, 8, 12];
        indices := p eq 2 select [1, 3, 4, 5] else [2, 4, 5];
        width := p eq 2 select 32 else 27;
        keys := [0, 1, p^width-1] cat [sample[1] : sample in model_samples(R)];
        assert #keys gt 0;
        for key in keys do
            W := aInvariantsFromKey(R, key);
            affine := [P | Evaluate(a, [z, P!1]) : a in W];
            infinity := [P | Evaluate(a, [P!1, z]) : a in W];
            for i in [1..5] do
                assert K!Evaluate(infinity[i], u) eq u^degrees[i]*Evaluate(affine[i], 1/u);
            end for;
            a1, a2, a3, a4, a6 := Explode([K | Evaluate(a, 1/u) : a in affine]);
            f := (y/u^6)^2+a1*(x/u^4)*(y/u^6)+a3*(y/u^6)
                -(x/u^4)^3-a2*(x/u^4)^2-a4*(x/u^4)-a6;
            b1, b2, b3, b4, b6 := Explode([K | Evaluate(a, u) : a in infinity]);
            g := y^2+b1*x*y+b3*y-x^3-b2*x^2-b4*x-b6;
            assert u^12*f eq g;
            assert FibrationKey(W) eq key;
        end for;
        for i in indices do
            for j in [0..degrees[i]] do
                W := [R | 0, 0, 0, 0, 0];
                W[i] := (p-1)*s^(degrees[i]-j)*t^j;
                key := (p-1)*p^(&+([0] cat [degrees[k]+1 : k in indices | k gt i])+degrees[i]-j);
                assert aInvariantsFromKey(R, key) eq W;
                assert Evaluate(W[i], [z, P!1]) eq (p-1)*z^(degrees[i]-j);
                assert Evaluate(W[i], [P!1, z]) eq (p-1)*z^j;
            end for;
        end for;
    end for;
end procedure;

// Encoding bounds are smaller than uint64; large integers must not wrap or truncate.
procedure test_model_key_validation()
    for p in [2, 3] do
        R := PolynomialRing(GF(p), 2);
        bound := p eq 2 select 2^32 else 3^27;
        assert bound lt 2^64;
        for key in [0, 1, bound-1] do
            assert FibrationKey(aInvariantsFromKey(R, key)) eq key;
        end for;
        for key in [-1, bound, bound+1, 2^64-1, 2^64, 2^64+1] do
            expect_error(procedure() W := aInvariantsFromKey(R, key); end procedure,
                "key is outside the encoding range");
        end for;
        expect_error(procedure() W := aInvariantsFromKey(R, "1"); end procedure,
            "Bad argument types");
        for W in [[R | 0, 0, 0], [R | R.1, 0, 0, 0, 0]] do
            expect_error(procedure() key := FibrationKey(W); end procedure,
                "five homogeneous coefficients");
            expect_error(procedure() kind, minimality := FibrationProperties(W); end procedure,
                "five homogeneous coefficients");
        end for;
        S := PolynomialRing(GF(p), 1);
        expect_error(procedure() W := aInvariantsFromKey(S, 0); end procedure,
            "Bad argument types");
        S := PolynomialRing(GF(p), 3);
        expect_error(procedure() W := aInvariantsFromKey(S, 0); end procedure,
            "binary polynomial ring");
    end for;
end procedure;

test_model_keys_and_cells();
test_model_charts();
test_model_key_validation();
Write(resultpath, "MODEL_MAGMA_OK");
quit;
