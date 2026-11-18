SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/support/_native_fixed_models.m": fixed_models;

function encode_fibres(rows)
    entries := [Strings() |];
    for r in rows do
        cycles := IsEmpty(r`ComponentCycles) select "[]"
            else Join([Sprint(n) : n in r`ComponentCycles], " ");
        fields := [Sprint(r`Place), Sprint(r`Degree), r`Type, Sprint(r`Index),
            Sprint(r`DiscriminantValuation), Sprint(r`Conductor),
            Sprint(r`Tamagawa), Sprint(r`Trace), cycles];
        Append(~entries, Join(fields, ":"));
    end for;
    return IsEmpty(entries) select "[]" else Join(entries, ";");
end function;

procedure check_known_aggregates(E, expected)
    result := FibrationLPolynomial(E);
    assert <result`Sign, result`Degree, result`Height,
        result`TwistHeight, result`TamagawaExponent> eq expected;
end procedure;

procedure test_known_local_aggregates()
    F<t> := FunctionField(GF(2));
    // Split I6, I6, I1: explicit discriminant and nodal fibres.
    check_known_aggregates(EllipticCurve([F | t, 0, t^6+1, 0, 0]),
        <-1, 5, 2, 4, 3>);
    // Hand Tate algorithm: I9* becomes III under Frobenius.
    check_known_aggregates(EllipticCurve([F | 0, 0, 1, t^6, 0]),
        <1, 7, 2, 1, -1>);
    // Hand Tate algorithm: III becomes I21* under Frobenius.
    check_known_aggregates(EllipticCurve([F | 0, 0, 1, t^7, 0]),
        <-1, 19, 2, 4, 1>);
    F<t> := FunctionField(GF(3));
    // Four degree-4 nonsplit nodes: explicit discriminant factorization.
    check_known_aggregates(EllipticCurve([F | 0, t^3+t, 0, t^8+1, 0]),
        <1, 12, 2, 6, 0>);
    // Hand Tate algorithm: wild II before and after Frobenius.
    check_known_aggregates(EllipticCurve([F | 0, 0, 0, t^8, t+1]),
        <1, 20, 2, 2, 0>);
end procedure;

procedure write_local_reference(path)
    empty := Open(path, "w");
    delete empty;
    fprintf path, "q,key,fibres,twist_fibres,degree,sign,height,twist_height,tamagawa_exponent\n";
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        fixed := fixed_models(R);
        models := fixed;
        // Full periods modulo 2^32 and 3^27 require multiplier 1 mod 12
        // and an increment coprime to 6.
        key := 17;
        sampled_keys := {Integers() |};
        for i in [1..100] do
            key := (1103515249*key+12347) mod p^(p eq 2 select 32 else 27);
            Include(~sampled_keys, key);
            W := aInvariantsFromKey(R, key);
            kind, minimality := FibrationProperties(W);
            if kind eq "smooth" and minimality eq "everywhere minimal" then
                Append(~models, W);
            end if;
        end for;
        assert #sampled_keys eq 100;
        for W in models do
            E := EllipticCurve([F | Evaluate(a, [t, F!1]) : a in W]);
            original := FibrationLocalData(E);
            twisted := FibrationLocalData(E : FrobeniusTwist := true);
            height := &+[r`Degree*r`DiscriminantValuation : r in original] div 12;
            twist_height := &+[r`Degree*r`DiscriminantValuation : r in twisted] div 12;
            exponent := &+[Valuation(r`Tamagawa, p) : r in twisted]
                - &+[Valuation(r`Tamagawa, p) : r in original];
            sign := IsEven(twist_height-height+exponent) select 1 else -1;
            degree := &+[r`Degree*r`Conductor : r in original]-4;
            for chart in [W, [R | Evaluate(a, [u, s]) : a in W]] do
                E := EllipticCurve([F | Evaluate(a, [t, F!1]) : a in chart]);
                original := FibrationLocalData(E);
                twisted := FibrationLocalData(E : FrobeniusTwist := true);
                fields := [Sprint(p), Sprint(FibrationKey(chart)), encode_fibres(original),
                    encode_fibres(twisted), Sprint(degree), Sprint(sign), Sprint(height),
                    Sprint(twist_height), Sprint(exponent)];
                fprintf path, "%o\n", Join(fields, ",");
            end for;
        end for;
    end for;
end procedure;

test_known_local_aggregates();
write_local_reference(referencepath);
quit;
