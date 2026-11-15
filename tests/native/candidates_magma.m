SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/_support.m": expect_error;

// The Python prototypes supply literal lists, including the Hermite exclusions.
procedure check_candidates(E, expected)
    data := FibrationInvariants(E);
    sha, regulators := FibrationCandidates(data, FibrationLocalData(E));
    assert sha eq expected;
    assert #sha eq #regulators and #sha gt 0;
    assert sha eq Sort(Setseq(Seqset(sha)));
    for i in [1..#sha] do
        assert sha[i] gt 0 and IsSquare(sha[i]);
        assert regulators[i] gt 0;
        assert sha[i]*regulators[i] eq data`ShaRegulatorProduct;
    end for;
end procedure;

// Actual K3 equations retain both rank-one and higher-rank ambiguity.
procedure test_fixed_candidates()
    F<t> := FunctionField(GF(2));
    check_candidates(EllipticCurve([F | 0, 0, 1, t^6, 0]), [1]);
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    check_candidates(E, [1, 4]);
    sha, regulators := FibrationCandidates(E);
    assert sha eq [1, 4] and regulators eq [Rationals() | 14, 7/2];
    check_candidates(EllipticCurve([F | t, 0, t^6+1, 0, 0]), [1]);
    E := EllipticCurve([F | t^2, 0, 1, 0, 0]);
    check_candidates(E, [1]);
    sha, regulators := FibrationCandidates(E);
    assert sha eq [1] and regulators eq [Rationals() | 1];
    F<t> := FunctionField(GF(3));
    check_candidates(EllipticCurve([F | 0, 0, 0, t^8, t+1]), [1]);
    check_candidates(EllipticCurve([F | 0, t^3+t, 0, t^8+1, 0]), [1, 4]);
    check_candidates(EllipticCurve([F | 0, 0, 0, 1, t^7]), [1, 4, 9]);
    E := EllipticCurve([F | 0, t^4+1, 0, t^4, 0]);
    check_candidates(E, [1]);
    sha, regulators := FibrationCandidates(E);
    assert regulators eq [Rationals() | 3/4];
end procedure;

// Supplied arithmetic controls pin every rational-component correction family.
procedure test_rational_component_corrections()
    D := recformat<Rank, PicardNumber, ArtinTateProduct, ShaRegulatorProduct>;
    F := recformat<Type, Index, Tamagawa, Degree>;
    cases := [
        <"I", 1, 1, 1, [Rationals() | 4]>,
        <"I", 4, 4, 1, [Rationals() | 4, 13/4, 3]>,
        <"I", 4, 2, 1, [Rationals() | 4, 3]>,
        <"I", 3, 1, 1, [Rationals() | 4]>,
        <"I*", 0, 1, 1, [Rationals() | 4]>,
        <"I*", 0, 2, 1, [Rationals() | 4, 3]>,
        <"I*", 0, 4, 1, [Rationals() | 4, 3]>,
        <"I*", 9, 2, 1, [Rationals() | 4, 3]>,
        <"I*", 9, 4, 1, [Rationals() | 4, 3, 11/4, 3/4]>,
        <"II", 0, 1, 1, [Rationals() | 4]>,
        <"II*", 0, 1, 1, [Rationals() | 4]>,
        <"III", 0, 2, 1, [Rationals() | 4, 7/2]>,
        <"III", 0, 2, 3, [Rationals() | 4, 5/2]>,
        <"IV", 0, 1, 1, [Rationals() | 4]>,
        <"IV", 0, 3, 1, [Rationals() | 4, 10/3]>,
        <"IV*", 0, 1, 1, [Rationals() | 4]>,
        <"IV*", 0, 3, 1, [Rationals() | 4, 8/3]>,
        <"III*", 0, 2, 1, [Rationals() | 4, 5/2]>];
    heights := [Rationals() | 4, 7/2, 10/3, 13/4, 3, 11/4, 8/3, 5/2, 3/4];
    for c in cases do
        rows := [rec<F | Type := c[1], Index := c[2], Tamagawa := c[3], Degree := c[4]>];
        for h in heights do
            data := rec<D | Rank := 1, PicardNumber := 2,
                ArtinTateProduct := 1, ShaRegulatorProduct := h>;
            sha, regulators := FibrationCandidates(data, rows);
            assert sha eq (h in c[5] select [Integers() | 1] else [Integers() |]);
            assert regulators eq (h in c[5] select [Rationals() | h] else [Rationals() |]);
        end for;
    end for;
end procedure;

// Large square divisors must be enumerated exactly without a square-root scan.
procedure test_large_exact_candidates()
    D := recformat<Rank, PicardNumber, ArtinTateProduct, ShaRegulatorProduct>;
    data := rec<D | Rank := 1, PicardNumber := 3,
        ArtinTateProduct := 4*9^20, ShaRegulatorProduct := 4*9^20>;
    sha, regulators := FibrationCandidates(data, []);
    assert sha eq [9^i : i in [0..20]];
    assert regulators eq [Rationals() | 4*9^(20-i) : i in [0..20]];
    assert sha[#sha] gt 2^53;
    assert &and[sha[i]*regulators[i] eq 4*9^20 : i in [1..#sha]];
    data`Rank := 0;
    data`PicardNumber := 2;
    data`ShaRegulatorProduct := 9^20;
    sha, regulators := FibrationCandidates(data, []);
    assert sha eq [9^20] and regulators eq [Rationals() | 1];
    data`ShaRegulatorProduct := 2;
    sha, regulators := FibrationCandidates(data, []);
    assert IsEmpty(sha) and IsEmpty(regulators);
end procedure;

// Above the tabulated Hermite ranks, the conservative bound retains valid pairs.
procedure test_high_rank_hermite_bound()
    D := recformat<Rank, PicardNumber, ArtinTateProduct, ShaRegulatorProduct>;
    for r in [3..20] do
        data := rec<D | Rank := r, PicardNumber := r+2,
            ArtinTateProduct := 4^r, ShaRegulatorProduct := 4^r>;
        sha, regulators := FibrationCandidates(data, []);
        assert 1 in sha;
        assert &and[sha[i]*regulators[i] eq 4^r : i in [1..#sha]];
        if r eq 3 then assert sha eq [1]; end if;
        if r in [4, 5, 6] then assert sha eq [1, 4]; end if;
    end for;
    data := rec<D | Rank := 21, PicardNumber := 22,
        ArtinTateProduct := 1, ShaRegulatorProduct := 1>;
    expect_error(procedure() sha, regulators := FibrationCandidates(data, []);
        end procedure, "rank must be between 0 and 20");
end procedure;

test_fixed_candidates();
test_rational_component_corrections();
test_large_exact_candidates();
test_high_rank_hermite_bound();
Write(resultpath, "CANDIDATES_MAGMA_OK");
quit;
