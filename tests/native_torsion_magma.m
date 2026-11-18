import "tests/support/_support.m": expect_error;

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

test_five_torsion_points();
test_torsion_domain();
