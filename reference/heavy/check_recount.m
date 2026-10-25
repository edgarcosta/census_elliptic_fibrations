SetQuitOnError(true);
AttachSpec("spec");
load "reference/f2/elliptic_fibration_orbits_f2.m";
load "reference/f2/k3_minimality_f2.m";
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/f3/ef3_stream.m";
load "reference/support/classification.m";
load "reference/support/f2_prefix_match.m";
load "reference/support/f3_prefix_match.m";
load "reference/support/heavy.m";
load "reference/support/parallel.m";

// Compares every cell on every prefix with independent scanned orbits.
procedure test_full_recount(shard,processes)
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        new := FibrationPrefixes(R);
        old := p eq 2 select FibrationOrbitPrefixesF2(R) else FibrationOrbitPrefixesF3(R);
        match := p eq 2 select reference_f2_prefix_match(R,old,new)
            else reference_f3_prefix_match(R,old,new);
        indices := heavy_prefix_indices(#new);
        selected := [ indices[i] : i in [1..#indices] | (i-1) mod processes eq shard-1 ];
        if IsEmpty(selected) then continue; end if;
        n,cells,rows := NumberOfFibrations(R : PrefixIndices := selected);
        assert [ row[1] : row in rows ] eq selected;
        for row in rows do
            j := row[1];
            expected := [ row[2][<kind,minimum>] : minimum in
                ["everywhere minimal","nonminimal somewhere"], kind in
                ["smooth","regular cusp","nonregular"] ];
            P := old[match[j]];
            actual := [0,0,0,0,0,0];
            points := 0;
            if p eq 2 or (P[1] eq 0 and P[2] eq 0) then
                if p eq 2 then
                    keys,info := FibrationFiberOrbitKeysF2(P : GenericSmoothOnly:=false);
                else
                    keys,info := FibrationFiberOrbitKeysF3(P : GenericSmoothOnly:=false);
                end if;
                points := p^info[1];
                assert #keys eq info[3];
                for key in keys do
                    W := p eq 2 select FibrationFromKeyF2(R,key) else EF3_AsLongModel(FibrationFromKeyF3(R,key));
                    if p eq 2 and FibrationDiscriminantF2(W) ne 0 then
                        c := IsEverywhereMinimalF2(W) select 1 else 2;
                    else
                        kind,minimum := reference_properties(W);
                        c := 2*(Index(["smooth","regular cusp","nonregular"],kind)-1)
                            +Index(["everywhere minimal","nonminimal somewhere"],minimum);
                    end if;
                    actual[c] +:= 1;
                end for;
            else
                stats := FibrationStreamPrefixF3(P,"unused" : WriteRepresentatives:=false);
                assert stats`Complete;
                Q,maps,S := EF3_FiberData(P);
                assert stats`PointsVisited eq 3^Dimension(Q);
                points := stats`PointsVisited;
                actual[1] := stats`K3Orbits;
                actual[2] := stats`NonMinimalOne+stats`NonMinimalSeveral;
                assert stats`SmoothOrbits eq actual[1]+actual[2];
                if P[1] ne 0 and P[2] eq 0 then
                    kind,minimum := reference_properties([R|0,P[1],0,0,0]);
                    assert kind eq "nonregular";
                    actual[minimum eq "everywhere minimal" select 5 else 6] := 1;
                end if;
                assert &+actual eq stats`AllOrbits;
            end if;
            for c in [1..6] do
                error if actual[c] ne expected[c], Sprintf(
                    "F%o public prefix %o reference prefix %o cell %o: package %o reference %o",
                    p,j,match[j],c,expected[c],actual[c]);
            end for;
        end for;
    end for;
end procedure;

try
    processes := heavy_processes();
    value := GetEnv("REFERENCE_SHARD");
    if value ne "" then
        test_full_recount(heavy_shard(processes),processes);
    elif processes eq 1 then
        test_full_recount(1,1);
    else
        reference_parallel_check("reference/heavy/check_recount.m",processes);
    end if;
catch e
    printf "MISMATCH: full recount %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
