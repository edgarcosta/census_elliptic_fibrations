SetQuitOnError(true);
AttachSpec("spec");
import "tests/_support.m": expect_error;
import "tests/_native_torsion_samples.m": torsion_samples;

// Exact equality with Magma pins completeness, orders, and the group law.
procedure check_torsion(E, expected)
    factors, points := FibrationTorsion(E);
    T, mp := TorsionSubgroup(E);
    assert factors eq expected and factors eq Invariants(T);
    assert #points eq #T and #Seqset(points) eq #points;
    assert IsId(points[1]);
    elements := [s : s in T];
    images := [mp(s) : s in elements];
    assert Seqset(points) eq Seqset(images);
    for P in points do
        assert Height(P) eq 0;
        s := elements[Index(images, P)];
        assert Order(P) eq Order(s);
        assert -P in points;
        for Q in points do
            assert P+Q in points;
            assert mp(s+elements[Index(images, Q)]) eq P+Q;
        end for;
    end for;
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

// A denominator of degree one and its reciprocal require the full bound m=1.
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
test_torsion_at_the_section_bound();
test_rank_zero_torsion();
test_torsion_domain();
Write(resultpath, "TORSION_MAGMA_OK");
quit;
