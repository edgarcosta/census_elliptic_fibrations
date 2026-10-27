SetQuitOnError(true);
AttachSpec("spec");
load "reference/support/heavy.m";
load "reference/support/parallel.m";

// Every written key is checked for distinctness, round trip and selected cell.
procedure test_prefix_listing(R,j,expected)
    p := Characteristic(R);
    kinds := ["smooth","regular cusp","nonregular"];
    minima := ["everywhere minimal","nonminimal somewhere"];
    scratch := GetCurrentDirectory() cat "/artifacts/test-tmp";
    failures := 0;
    for c in [1..6] do
        kind := kinds[(c+1) div 2];
        minimum := minima[2-(c mod 2)];
        stem,handle := TemporaryFile("listing-" : ParentDirectory := scratch);
        delete handle; Remove(stem);
        binary_path := Sprintf("%o_%o.bin",stem,IntegerToString(10000+j)[2..5]);
        stats_path := Sprintf("%o_%o.stats",stem,IntegerToString(10000+j)[2..5]);
        try
            n := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices:=[j],
                GenericFibreKinds:=[kind],Minimalities:=[minimum]);
            if n ne expected[c] then
                printf "MISMATCH,F%o,prefix=%o,cell=%o,count=%o,written=%o\n",p,j,c,expected[c],n;
                failures +:= 1;
            end if;
            F := Open(binary_path,"rb");
            seen := {Integers()|};
            checked := 0;
            read_count := 0;
            while true do
                bytes := ReadBytes(F,32768);
                if #bytes eq 0 then break; end if;
                assert #bytes mod 8 eq 0;
                for offset in [0..(#bytes div 8)-1] do
                    key := &+[bytes[8*offset+b+1]*256^b:b in [0..7]];
                    read_count +:= 1;
                    if key in seen then
                        printf "MISMATCH,F%o,prefix=%o,cell=%o,duplicate_key=%o\n",p,j,c,key;
                        failures +:= 1;
                    end if;
                    Include(~seen,key);
                    W := aInvariantsFromKey(R,key);
                    back := FibrationKey(W);
                    actual_kind,actual_minimum := FibrationProperties(W);
                    if back ne key or actual_kind ne kind or actual_minimum ne minimum then
                        printf "MISMATCH,F%o,prefix=%o,cell=%o,key=%o,roundtrip=%o,expected=%o|%o,actual=%o|%o\n",
                            p,j,c,key,back,kind,minimum,actual_kind,actual_minimum;
                        failures +:= 1;
                    end if;
                    checked +:= 1;
                end for;
            end while;
            delete F;
            assert read_count eq n and #seen eq n and checked eq n;
            assert "complete=1\n" in Read(stats_path);
            Remove(binary_path);
            Remove(stats_path);
            Remove(stem cat ".csv");
        catch e
            printf "MISMATCH,F%o,prefix=%o,cell=%o\n%o\n",p,j,c,e;
            failures +:= 1;
            for path in [binary_path,stats_path,stem cat ".csv"] do
                if FileExists(path) then Remove(path); end if;
            end for;
        end try;
    end for;
    error if failures ne 0, "MISMATCH: representative listing failed";
end procedure;

procedure test_full_listing(shard,processes)
    assert System("mkdir -p artifacts/test-tmp") eq 0;
    for p in [2,3] do
        R := PolynomialRing(GF(p),2);
        indices := heavy_prefix_indices(#FibrationPrefixes(R));
        selected := [ indices[i] : i in [1..#indices] | (i-1) mod processes eq shard-1 ];
        if IsEmpty(selected) then continue; end if;
        n,cells,rows := NumberOfFibrations(R : PrefixIndices := selected);
        assert [ row[1] : row in rows ] eq selected;
        for row in rows do
            expected := [ row[2][<kind,minimum>] : minimum in
                ["everywhere minimal","nonminimal somewhere"], kind in
                ["smooth","regular cusp","nonregular"] ];
            test_prefix_listing(R,row[1],expected);
        end for;
    end for;
end procedure;

try
    processes := heavy_processes();
    value := GetEnv("REFERENCE_SHARD");
    if value ne "" then
        test_full_listing(heavy_shard(processes),processes);
    elif processes eq 1 then
        test_full_listing(1,1);
    else
        reference_parallel_check("reference/heavy/check_listing.m",processes);
    end if;
catch e
    printf "MISMATCH: full listing %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
