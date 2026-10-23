// Complete small-fiber comparisons, binary/text round trips, and bounded
// C!=0 smoke tests. This is not a census benchmark.
load "ef_common.m";
load "ef3_models.m";
load "ef3_orbits.m";
load "k3_minimality_f3.m";
load "ef3_stream.m";

procedure EF3_RunStreamChecks()
R<s,t> := PolynomialRing(GF(3),2);
P := FibrationOrbitPrefixesF3(R);
tested := {Integers()|};
for prefix in P do
    if prefix[1] ne 0 or prefix[2] eq 0 then continue; end if;
    Q,maps,S := EF3_FiberData(prefix);
    d := Dimension(Q);
    if d in tested then continue; end if;
    reference,info := FibrationFiberOrbitKeysF3(prefix);
    expected := [Integers()|];
    for key in reference do
        ok,bad := IsEverywhereMinimalF3(FibrationFromKeyF3(R,key));
        if ok then Append(~expected,key); end if;
    end for;
    filename := Sprintf("stream_check_%o.bin",d);
    stats := FibrationStreamPrefixF3(prefix,filename);
    assert stats`Complete and stats`K3Orbits eq #expected;
    F := Open(filename,"rb");
    actual := [Integers()|];
    while true do
        found,key := EF3_ReadBinaryKey(F);
        if not found then break; end if;
        Append(~actual,key);
    end while;
    delete F;
    assert actual eq expected;
    filename := Sprintf("stream_check_%o.keys",d);
    stats := FibrationStreamPrefixF3(prefix,filename : BinaryOutput:=false);
    F := Open(filename,"r");
    actual := [Integers()|];
    while true do
        line := Gets(F);
        if IsEof(line) then break; end if;
        Append(~actual,StringToInteger(line));
    end while;
    delete F;
    assert actual eq expected;
    printf "Complete d=%o fiber: %o K3 representatives; binary/text match\n",d,#actual;
    Include(~tested,d);
end for;
assert tested eq {8,9};

// Exercise both trivial/nontrivial maps and auto/non-auto minimality.
testedKinds := {};
for prefix in P do
    if prefix[1] eq 0 then continue; end if;
    Q,maps,S := EF3_FiberData(prefix);
    automatic := K3F3_PrefixAutomaticallyMinimal(prefix[1],prefix[2]);
    kind := <#maps eq 1,automatic>;
    if kind in testedKinds then continue; end if;
    filename := Sprintf("stream_check_C_%o.bin",#testedKinds);
    stats := FibrationStreamPrefixF3(prefix,filename : PointLimit:=2000);
    assert not stats`Complete and stats`PointsVisited eq 2000;
    F := Open(filename,"rb");
    seen := {Integers()|};
    while true do
        found,key := EF3_ReadBinaryKey(F);
        if not found then break; end if;
        assert key notin seen;
        Include(~seen,key);
        model := FibrationFromKeyF3(R,key);
        assert model[1] eq prefix[1] and model[2] eq prefix[2];
        ok,bad := IsEverywhereMinimalF3(model);
        assert ok;
    end while;
    delete F;
    assert #seen eq stats`K3Orbits;
    printf "C!=0 bounded check <trivial,auto>=%o: %o K3 keys\n",kind,#seen;
    Include(~testedKinds,kind);
    if #testedKinds eq 4 then break; end if;
end for;
assert #testedKinds eq 4;

// Test summary and completion-sidecar writes on a small j=0 prefix.
i := Min([j:j in [1..#P] | P[j][1] eq 0 and P[j][2] ne 0]);
result := CensusF3(R,[i],"stream_check_driver");
assert #result eq 1 and result[1][2]`Complete;
assert i eq 2;
present,F := OpenTest("stream_check_driver_0002.stats","r");
assert present;
delete F;
print "Streaming checks passed.";
end procedure;

EF3_RunStreamChecks();
quit;
