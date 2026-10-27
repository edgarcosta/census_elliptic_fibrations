AttachSpec("spec");
load "reference/f2/elliptic_fibration_orbits_f2.m";
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/support/classification.m";
load "reference/support/singular.m";

// Recounts every singular orbit without the package's orbit engine or loci.
procedure test_singular_classification()
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        models,cusp_count := reference_singular_models(R);
        cells := reference_singular_cells(models);
        expected := p eq 2 select [680,38,7,11] else [294,20,12,12];
        assert cells eq expected;
        assert #models eq (p eq 2 select 736 else 338);
        n,table := NumberOfFibrations(R : GenericFibreKinds := ["regular cusp","nonregular"]);
        assert n eq #models;
        assert [ table[<kind,minimum>] : minimum in ["everywhere minimal","nonminimal somewhere"],
            kind in ["regular cusp","nonregular"] ] eq cells;
        overlaps := 0;
        stabilizer_checks := 0;
        // Reference transformations cover the full prefix stabilizers.
        // Sextic shifts over F2 and cube shifts over F3 are already quotiented out.
        for W in models do
            check_reference_properties(W);
            kind,minimum,certificate := reference_properties(W);
            if certificate[1] eq "finite" and #certificate[2] gt 1 then
                overlaps +:= 1;
            end if;
            if p eq 2 then
                for g in GL(2,GF(2)), n in [0..7] do
                    B := TransformNormalizedF2(W,Matrix(g),EF2FormFromMask(R,2,n),R!0);
                    if B[1] ne W[1] or B[3] ne W[3] then continue; end if;
                    next_kind,next_minimum := FibrationProperties(B);
                    assert <next_kind,next_minimum> eq <kind,minimum>;
                    stabilizer_checks +:= 1;
                end for;
            else
                for M in EF3_BaseMatrices(GF(3)) do
                    C := EF_BaseChange(W[2],M);
                    A := EF_BaseChange(W[4],M);
                    if C ne W[2] then continue; end if;
                    if C eq 0 then
                        if A ne W[4] then continue; end if;
                        q := R!0;
                    else
                        divisible,q := IsDivisibleBy(A-W[4],C);
                        if not divisible then continue; end if;
                    end if;
                    B := EF3_AsLongModel(EF3_Transform([W[2],W[4],W[5]],M,q));
                    next_kind,next_minimum := FibrationProperties(B);
                    assert <next_kind,next_minimum> eq <kind,minimum>;
                    stabilizer_checks +:= 1;
                end for;
            end if;
        end for;
        assert overlaps gt 0 and stabilizer_checks gt #models;
        printf "F%o singular: regular cusp minimal=%o nonminimal=%o; nonregular minimal=%o nonminimal=%o; total=%o (cuspidal=%o); overlaps=%o; stabilizer checks=%o\n",
            p,cells[1],cells[2],cells[3],cells[4],#models,cusp_count,overlaps,stabilizer_checks;
    end for;
end procedure;

test_singular_classification();
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
