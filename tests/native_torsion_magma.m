import "tests/support/_support.m": expect_error;
import "tests/support/_native_torsion_samples.m": torsion_samples;

procedure check_torsion(E, expected)
    factors, points := FibrationTorsion(E);
    assert factors eq expected;
    assert #Seqset(points) eq #points;
    assert IsId(points[1]);
end procedure;

// Both characteristics include trivial, cyclic, and characteristic-primary torsion.
procedure test_fixed_torsion_groups()
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        for sample in torsion_samples(F) do
            equation, expected := Explode(sample);
            check_torsion(EllipticCurve(equation), expected);
        end for;
    end for;
end procedure;

// Explicit multiples of (0,t^2), good reduction and the 3-division polynomial give C5.
procedure test_five_torsion_points()
    F<t> := FunctionField(GF(3));
    a := [F | 0, 1+t^4, 0, t^2-t^4, t^4];
    coordinates := [[F | 0, t^2], [F | t^2, 2*t^4],
        [F | t^2, t^4], [F | 0, 2*t^2]];
    weights := [2, 4, 6, 8, 12];
    for reciprocal in [false, true] do
        equation := reciprocal select [F | t^weights[i]*Evaluate(a[i], 1/t)
            : i in [1..5]] else a;
        E := EllipticCurve(equation);
        factors, points := FibrationTorsion(E);
        expected := [E!0] cat [E!(reciprocal select
            [t^4*Evaluate(P[1], 1/t), t^6*Evaluate(P[2], 1/t)] else P)
            : P in coordinates];
        assert factors eq [5] and Seqset(points) eq Seqset(expected);
        assert #points eq 5;
        assert [Order(P) : P in expected] eq [1, 5, 5, 5, 5];
    end for;
end procedure;

procedure test_torsion_at_the_section_bound()
    F<t> := FunctionField(GF(2));
    samples := torsion_samples(F);
    E := EllipticCurve(samples[4][1]);
    factors, points := FibrationTorsion(E);
    assert E![1/t^2, 1/t^3] in points;
    E := EllipticCurve(samples[5][1]);
    factors, points := FibrationTorsion(E);
    assert E![t^6, t^9] in points;
end procedure;

// Rank zero includes both trivial torsion and nonidentity torsion points.
procedure test_rank_zero_torsion()
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        samples := torsion_samples(F);
        E := EllipticCurve(samples[#samples][1]);
        assert Evaluate(LFunction(E), 1/p) ne 0;
    end for;
    F<t> := FunctionField(GF(2));
    E := EllipticCurve(torsion_samples(F)[6][1]);
    assert Evaluate(LFunction(E), 1/2) ne 0;
end procedure;

// The K3 proof requires the supported field and an everywhere minimal chart.
procedure test_torsion_domain()
    F<t> := FunctionField(GF(4));
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    expect_error(procedure() factors := FibrationTorsion(E); end procedure,
        "constant field must be GF(2) or GF(3)");
    F<t> := FunctionField(GF(2));
    for a in [[F | 0, 0, 1, 0, 0], [F | 0, 0, t^6, t^8, 0]] do
        E := EllipticCurve(a);
        expect_error(procedure() factors := FibrationTorsion(E); end procedure,
            "equation must define an everywhere minimal elliptic K3");
    end for;
    for a in [[F | 0, 0, 1, t^9, 0], [F | 0, 0, 1, 1/t, 0]] do
        E := EllipticCurve(a);
        expect_error(procedure() factors := FibrationTorsion(E); end procedure,
            "coefficients must be polynomials");
    end for;
end procedure;

test_fixed_torsion_groups();
test_five_torsion_points();
test_torsion_at_the_section_bound();
test_rank_zero_torsion();
test_torsion_domain();
