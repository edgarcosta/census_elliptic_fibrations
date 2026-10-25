SetQuitOnError(true);
AttachSpec("spec");
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/f3/ef3_stream.m";
load "reference/support/f3_prefix_match.m";

// Recounts the nonautomatic C!=0 fibres omitted by check_gmodule.m.
procedure test_f3_remaining_smooth_splits(shard, processes)
    R<s,t> := PolynomialRing(GF(3),2);
    old := FibrationOrbitPrefixesF3(R);
    new := FibrationPrefixes(R);
    assert #old eq 1860;
    match := reference_f3_prefix_match(R,old,new);
    orders := {};
    indices := [];
    for j in [1..#new] do
        P := old[match[j]];
        if P[1] eq 0 or K3F3_PrefixAutomaticallyMinimal(P[1],P[2]) then
            continue;
        end if;
        Q,maps,S := EF3_FiberData(P);
        if #maps gt 1 and #maps notin orders then
            Include(~orders,#maps);
        else
            Append(~indices,j);
        end if;
    end for;
    assert #orders eq 4 and #indices eq 39;
    selected := [ indices[i] : i in [1..#indices] | (i-1) mod processes eq shard-1 ];
    assert not IsEmpty(selected);
    n,cells,rows := NumberOfFibrations(R : PrefixIndices := selected);
    for row in rows do
        j := row[1]; c := row[2];
        stats := FibrationStreamPrefixF3(old[match[j]],"unused" :
            WriteRepresentatives := false);
        assert stats`Complete and stats`PointsVisited eq 3^13;
        error if stats`AllOrbits ne &+[ c[cell] : cell in Keys(c) ],
            Sprintf("F3 prefix %o: all-orbit count differs",j);
        smooth_minimal := c[<"smooth","everywhere minimal">];
        smooth_nonminimal := c[<"smooth","nonminimal somewhere">];
        error if stats`SmoothOrbits ne smooth_minimal+smooth_nonminimal,
            Sprintf("F3 prefix %o: smooth count differs",j);
        error if stats`K3Orbits ne smooth_minimal,
            Sprintf("F3 prefix %o: smooth minimal count differs",j);
        error if stats`NonMinimalOne+stats`NonMinimalSeveral ne smooth_nonminimal,
            Sprintf("F3 prefix %o: smooth nonminimal count differs",j);
    end for;
end procedure;

load "reference/support/parallel.m";

value := GetEnv("REFERENCE_PROCESSES");
processes := value eq "" select 8 else StringToInteger(value);
error if processes lt 1, "REFERENCE_PROCESSES must be positive";
processes := Min(processes,39);
value := GetEnv("REFERENCE_SHARD");
if value ne "" then
    shard := StringToInteger(value);
    assert shard ge 1 and shard le processes;
    test_f3_remaining_smooth_splits(shard,processes);
elif processes eq 1 then
    test_f3_remaining_smooth_splits(1,1);
else
    reference_parallel_check("reference/checks/check_smooth_f3.m",processes);
end if;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
