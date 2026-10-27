SetQuitOnError(true);
AttachSpec("spec");
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/support/f3_prefix_match.m";

// Homogeneous factors check the shortcut at finite places and infinity.
procedure test_f3_automatic_minimality()
    R<s,t> := PolynomialRing(GF(3),2);
    old := FibrationOrbitPrefixesF3(R);
    new := FibrationPrefixes(R);
    assert #old eq 1860;
    match := reference_f3_prefix_match(R,old,new);
    n,cells,rows := NumberOfFibrations(R);
    automatic_count := 0;
    smooth_total := 0;
    nodes := 0;
    for j in [1..#new] do
        P := old[match[j]];
        C := P[1]; A := P[2];
        if C eq 0 then continue; end if;
        can_be_bad := false;
        for factor in Factorization(C) do
            pi := factor[1]; m := factor[2];
            if m ge 2 and IsDivisibleBy(A,pi^Min(m,4)) then
                can_be_bad := true;
            end if;
        end for;
        assert K3F3_PrefixAutomaticallyMinimal(C,A) eq not can_be_bad;
        if can_be_bad then continue; end if;
        automatic_count +:= 1;
        total := EF3_FiberBurnsideCount(P);
        singular := A eq 0 select 1 else 0;
        c := rows[j][2];
        assert c[<"smooth","everywhere minimal">] eq total-singular;
        assert c[<"smooth","nonminimal somewhere">] eq 0;
        assert c[<"nonregular","everywhere minimal">] eq singular;
        assert c[<"nonregular","nonminimal somewhere">] eq 0;
        assert c[<"regular cusp","everywhere minimal">] eq 0;
        assert c[<"regular cusp","nonminimal somewhere">] eq 0;
        smooth_total +:= total-singular;
        nodes +:= singular;
    end for;
    assert automatic_count eq 890;
    assert smooth_total eq 1254370308 and nodes eq 12;
end procedure;

test_f3_automatic_minimality();
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
