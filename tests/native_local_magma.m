import "tests/support/_support.m": expect_error;

function local_tuple(r)
    return <r`Place, r`Degree, r`Type, r`Index, r`DiscriminantValuation,
        r`Conductor, r`Tamagawa, r`Trace, r`ComponentCycles>;
end function;

function component_count(r)
    case r`Type:
        when "I": return r`Index;
        when "I*": return r`Index+5;
        when "II": return 1;
        when "III": return 2;
        when "IV": return 3;
        when "IV*": return 7;
        when "III*": return 8;
        when "II*": return 9;
    end case;
    error "unexpected Kodaira family";
end function;

// Ogg's identity and component degrees apply even to wild or non-K3 twists.
procedure check_identities(rows)
    for r in rows do
        m := component_count(r);
        assert r`Conductor eq r`DiscriminantValuation+1-m;
        assert &+([0] cat r`ComponentCycles) eq r`Degree*(m-1);
        assert r`ComponentCycles eq Sort(r`ComponentCycles);
        assert &and[n mod r`Degree eq 0 : n in r`ComponentCycles];
    end for;
end procedure;

// Reciprocal charts must give identical Tate data at corresponding places.
procedure check_charts(E)
    F := BaseRing(E);
    t := F.1;
    p := Characteristic(F);
    a := aInvariants(E);
    weights := [2, 4, 6, 8, 12];
    reversed := EllipticCurve([F | t^weights[i]*Evaluate(a[i], 1/t) : i in [1..5]]);
    for twist in [false, true] do
        rows := FibrationLocalData(E : FrobeniusTwist := twist);
        other := FibrationLocalData(reversed : FrobeniusTwist := twist);
        assert #rows gt 0 and #rows eq #other;
        check_identities(rows);
        check_identities(other);
        P<z> := PolynomialRing(GF(p));
        for r in rows do
            if r`Place cmpeq "inf" then
                place := p;
            elif r`Place eq p then
                place := "inf";
            else
                f := P!IntegerToSequence(r`Place, p);
                f := P!Reverse(Coefficients(f));
                f /:= LeadingCoefficient(f);
                place := SequenceToInteger([Integers()!c : c in Coefficients(f)], p);
            end if;
            matched := [v : v in other | v`Place cmpeq place];
            assert #matched eq 1;
            v := matched[1];
            v`Place := r`Place;
            assert local_tuple(r) eq local_tuple(v);
        end for;
    end for;
end procedure;

// Literal nodal controls pin residue-degree cycles, traces and Tamagawa numbers.
procedure test_split_and_nonsplit_fibres()
    F<t> := FunctionField(GF(2));
    E := EllipticCurve([F | t, 0, t^6+1, 0, 0]);
    rows := FibrationLocalData(E);
    assert [local_tuple(r) : r in rows] eq [
        <3, 1, "I", 6, 6, 1, 6, 1, [1, 1, 1, 1, 1]>,
        <7, 2, "I", 6, 6, 1, 6, 1, [2, 2, 2, 2, 2]>,
        <73, 6, "I", 1, 1, 1, 1, 1, []>];
    check_charts(E);
    E := EllipticCurve([F | t, 0, 1, 0, t^12]);
    rows := FibrationLocalData(E);
    assert local_tuple(rows[1]) eq <262153, 18, "I", 1, 1, 1, 1, -1, []>;
    assert local_tuple(rows[2]) eq <"inf", 1, "III", 0, 6, 5, 2, 0, [1]>;
    check_charts(E);
    F<t> := FunctionField(GF(3));
    E := EllipticCurve([F | 0, t^3+t, 0, t^8+1, 0]);
    rows := FibrationLocalData(E);
    assert #rows eq 4;
    assert {<r`Degree, r`Index, r`Trace, r`Tamagawa> : r in rows}
        eq {<4, 1, -1, 1>, <4, 2, -1, 2>};
    twisted := FibrationLocalData(E : FrobeniusTwist := true);
    assert {<r`Index, r`ComponentCycles> : r in twisted}
        eq {<3, [8]>, <6, [4, 8, 8]>};
    check_charts(E);
    E := EllipticCurve([F | 0, t^4, 0, 1, t^12]);
    rows := FibrationLocalData(E);
    assert #rows eq 2 and &and[r`Trace eq 1 and r`Degree eq 12 : r in rows];
    check_charts(E);
end procedure;

// Tate residue quadratics and discriminant valuations distinguish split I* branches.
procedure test_positive_index_star_fibres()
    F<t> := FunctionField(GF(3));
    for sample in [<t^4, 1, 4, [1, 1, 1, 1, 1]>,
        <2*t^4, 1, 2, [1, 1, 1, 2]>, <2*t^5, 2, 4, [1 : i in [1..6]]>] do
        a6, n, c, cycles := Explode(sample);
        E := EllipticCurve([F | 0, t, 0, t^8, a6]);
        rows := [r : r in FibrationLocalData(E) | r`Place cmpeq 3];
        assert #rows eq 1;
        assert local_tuple(rows[1]) eq <3, 1, "I*", n, n+6, 2, c, 0, cycles>;
        rows := [r : r in FibrationLocalData(E : FrobeniusTwist := true)
            | r`Place cmpeq 3];
        assert #rows eq 1;
        twisted_cycles := c eq 4 select [1 : i in [1..3*n+4]]
            else [1 : i in [1..3*n+2]] cat [2];
        assert local_tuple(rows[1]) eq <3, 1, "I*", 3*n, 3*n+6, 2, c, 0, twisted_cycles>;
        check_charts(E);
    end for;
end procedure;

// The saved infinity controls require minimization after raising a_i to p.
procedure test_wild_frobenius_twists()
    F<t> := FunctionField(GF(2));
    for d in [6, 7] do
        E := EllipticCurve([F | 0, 0, 1, t^d, 0]);
        rows := FibrationLocalData(E);
        assert #rows eq 1;
        if d eq 6 then
            assert local_tuple(rows[1]) eq <"inf", 1, "I*", 9, 24, 11, 4, 0, [1 : i in [1..13]]>;
        else
            assert local_tuple(rows[1]) eq <"inf", 1, "III", 0, 24, 23, 2, 0, [1]>;
        end if;
        twisted := FibrationLocalData(E : FrobeniusTwist := true);
        assert #twisted eq 1;
        if d eq 6 then
            assert local_tuple(twisted[1]) eq <"inf", 1, "III", 0, 12, 11, 2, 0, [1]>;
        else
            assert local_tuple(twisted[1]) eq <"inf", 1, "I*", 21, 48, 23, 4, 0, [1 : i in [1..25]]>;
        end if;
        check_charts(E);
    end for;
    F<t> := FunctionField(GF(3));
    E := EllipticCurve([F | 0, 0, 0, t^8, t+1]);
    expected := <3, 1, "II", 0, 24, 24, 1, 0, []>;
    assert [local_tuple(r) : r in FibrationLocalData(E)] eq [expected];
    assert [local_tuple(r) : r in FibrationLocalData(E : FrobeniusTwist := true)] eq [expected];
    check_charts(E);
end procedure;

// Smooth minimal K3 models have discriminant divisor of degree 24.
procedure test_k3_discriminant_degrees()
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        forms := p eq 2 select [
            [R | s*u, 0, s^6+u^6, 0, 0],
            [R | 0, 0, u^6, s^6*u^2, 0],
            [R | 0, 0, u^6, s^7*u, 0]] else [
            [R | 0, s^3*u+s*u^3, 0, s^8+u^8, 0],
            [R | 0, 0, 0, s^8, (s+u)*u^11]];
        assert #forms gt 0;
        for W in forms do
            kind, minimality := FibrationProperties(W);
            assert kind eq "smooth" and minimality eq "everywhere minimal";
            E := EllipticCurve([F | Evaluate(a, [t, F!1]) : a in W]);
            rows := FibrationLocalData(E);
            assert &+[r`Degree*r`DiscriminantValuation : r in rows] eq 24;
            check_identities(rows);
        end for;
    end for;
end procedure;

// Constant good curves give no bad-place rows; unsupported fields are rejected.
procedure test_local_domain()
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        E := EllipticCurve(p eq 2 select [F | 0, 0, 1, 0, 0] else [F | 0, 0, 0, 1, 0]);
        assert IsEmpty(FibrationLocalData(E));
        assert IsEmpty(FibrationLocalData(E : FrobeniusTwist := true));
    end for;
    F<t> := FunctionField(GF(3));
    E := EllipticCurve([F | 0, t, 0, t^8, t^4]);
    // Scaling x by t^2 and y by t^3 recovers the polynomial equation.
    rational := EllipticCurve([F | 0, 1/t, 0, t^4, 1/t^2]);
    for twist in [false, true] do
        rows := FibrationLocalData(rational : FrobeniusTwist := twist);
        expected := FibrationLocalData(E : FrobeniusTwist := twist);
        assert [local_tuple(r) : r in rows] eq [local_tuple(r) : r in expected];
        at_zero := [r : r in rows | r`Place cmpeq 3];
        assert #at_zero eq 1;
        n := twist select 3 else 1;
        assert local_tuple(at_zero[1]) eq <3, 1, "I*", n, n+6, 2, 4, 0,
            [1 : i in [1..n+4]]>;
    end for;
    F<t> := FunctionField(GF(4));
    E := EllipticCurve([F | 0, 0, 1, t, 0]);
    expect_error(procedure() rows := FibrationLocalData(E); end procedure,
        "constant field must be GF(2) or GF(3)");
end procedure;

test_split_and_nonsplit_fibres();
test_positive_index_star_fibres();
test_wild_frobenius_twists();
test_k3_discriminant_degrees();
test_local_domain();
