// Streaming checks without file I/O.
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/f3/ef3_stream.m";

procedure EF3_RunStreamCoreChecks()
    R<s,t> := PolynomialRing(GF(3),2);
    P := FibrationOrbitPrefixesF3(R);
    tested := {Integers()|};
    for prefix in P do
        if prefix[1] ne 0 or prefix[2] eq 0 then continue; end if;
        Q,maps,S := EF3_FiberData(prefix);
        d := Dimension(Q);
        if d in tested then continue; end if;
        keys,info := FibrationFiberOrbitKeysF3(prefix);
        expected := 0; nonminimal := 0;
        for key in keys do
            ok,bad := IsEverywhereMinimalF3(FibrationFromKeyF3(R,key));
            if ok then expected +:= 1; else nonminimal +:= 1; end if;
            bytes := [(key div 256^j) mod 256:j in [0..7]];
            assert (&+[bytes[j+1]*256^j:j in [0..7]]) eq key;
        end for;
        stats := FibrationStreamPrefixF3(prefix,"unused" : WriteRepresentatives:=false);
        assert stats`Complete and stats`K3Orbits eq expected;
        assert stats`SmoothOrbits eq #keys;
        assert stats`NonMinimalOne+stats`NonMinimalSeveral eq nonminimal;
        printf "Complete d=%o stream agrees with reference: %o K3\n",d,expected;
        Include(~tested,d);
    end for;
    assert tested eq {8,9};
    testedKinds := {};
    for prefix in P do
        if prefix[1] eq 0 then continue; end if;
        Q,maps,S := EF3_FiberData(prefix);
        automatic := K3F3_PrefixAutomaticallyMinimal(prefix[1],prefix[2]);
        kind := <#maps eq 1,automatic>;
        if kind in testedKinds then continue; end if;
        stats := FibrationStreamPrefixF3(prefix,"unused" :
            WriteRepresentatives:=false,PointLimit:=2000);
        assert not stats`Complete and stats`PointsVisited eq 2000;
        assert stats`AutomaticMinimal eq automatic;
        assert stats`AllOrbits ge stats`SmoothOrbits;
        printf "C!=0 bounded stream <trivial,auto>=%o: %o K3\n",kind,stats`K3Orbits;
        Include(~testedKinds,kind);
        if #testedKinds eq 4 then break; end if;
    end for;
    assert #testedKinds eq 4;
    buckets := FibrationBenchmarkBucketsF3(P);
    autoCount := 0; otherCount := 0;
    autoPoints := 0; otherPoints := 0;
    for key in Keys(buckets) do
        if key[4] then
            autoCount +:= #buckets[key];
            autoPoints +:= #buckets[key]*3^key[2];
        else
            otherCount +:= #buckets[key];
            otherPoints +:= #buckets[key]*3^key[2];
        end if;
    end for;
    assert autoCount eq 1758 and otherCount eq 101;
    assert autoPoints eq 1424747394 and otherPoints eq 69015159;
    printf "Automatic/nonautomatic prefixes: %o / %o\n",autoCount,otherCount;
    printf "Automatic/nonautomatic scan points: %o / %o\n",autoPoints,otherPoints;
    print "Stream core checks passed; file I/O is NOT tested here.";
end procedure;

EF3_RunStreamCoreChecks();
Write(resultpath,"REFERENCE_CHECK_OK");
quit 0;
