// Streaming census. Load ef_common.m, ef3_models.m, ef3_orbits.m and
// k3_minimality_f3.m first.
// Output binary records are eight-byte unsigned little-endian model keys.
// No file header is written. The accompanying .stats/.csv files identify
// the prefix and counts. Text mode writes one decimal key per line.

EF3_RunStats := recformat<Dimension,MapCount,PointsVisited,AllOrbits,
    SmoothOrbits,K3Orbits,NonMinimalOne,NonMinimalSeveral,AutomaticMinimal,
    Complete,CPUSeconds>;

procedure EF3_RequireNewFile(filename)
    present,F := OpenTest(filename,"rb");
    if present then
        delete F;
        error "Output file already exists; use a new output stem: " cat filename;
    end if;
end procedure;

function EF3_ReadBinaryKey(F)
    bytes := ReadBytes(F,8);
    if #bytes eq 0 then return false,0; end if;
    assert #bytes eq 8;
    n := &+[ Integers() | bytes[i+1]*256^i : i in [0..7] ];
    assert n lt 3^27;
    return true,n;
end function;

function FibrationStreamPrefixF3(prefix, filename :
    WriteRepresentatives := true, BinaryOutput := true, PointLimit := 0,
    VerifyBurnside := true, Context := K3F3_default_context)
    // Returns a statistics record, never a sequence of representative keys.
    // Precondition: prefix was returned by FibrationOrbitPrefixesF3.
    // PointLimit=0 scans the complete fiber; a positive limit is for smoke
    // tests only, not for extrapolating census time from an initial segment.
    C := prefix[1]; A := prefix[2]; R := Parent(C); k := BaseRing(R);
    start := Cputime();
    Q,maps,S := EF3_FiberData(prefix);
    d := Dimension(Q);
    m := #maps;
    count := 3^d;
    if PointLimit gt 0 then count := Min(count,PointLimit); end if;
    complete := count eq 3^d;
    stats := rec< EF3_RunStats | Dimension:=d,MapCount:=m,PointsVisited:=0,
        AllOrbits:=0,SmoothOrbits:=0,K3Orbits:=0,NonMinimalOne:=0,
        NonMinimalSeveral:=0,AutomaticMinimal:=false,Complete:=complete,
        CPUSeconds:=0.0 >;
    if WriteRepresentatives then
        EF3_RequireNewFile(filename);
        F := Open(filename,BinaryOutput select "wb" else "w");
    end if;
    bytes := [ Integers() | ];
    lines := [ Strings() | ];
    buffered := 0;

    // The whole C=A=0 stratum is singular. Count it without enumeration.
    if C eq 0 and A eq 0 then
        stats`AllOrbits := EF3_AffineBurnsideCount(maps);
        stats`Complete := true;
        stats`CPUSeconds := Cputime(start);
        if WriteRepresentatives then Flush(F); delete F; end if;
        return stats;
    end if;

    automatic := K3F3_PrefixAutomaticallyMinimal(C,A : Context:=Context);
    stats`AutomaticMinimal := automatic;
    W := VectorSpace(k,d);
    lift := BasisMatrix(Q);
    if C ne 0 then
        assert d eq 13 and lift eq IdentityMatrix(k,13);
    end if;
    prefixPart := EF3_Key(EF_Coefficients(C,4) cat EF_Coefficients(A,8))*3^13;
    if m gt 1 then seen := [ false : i in [1..3^d] ]; end if;
    if m eq 1 then
        assert maps[1][1] eq IdentityMatrix(k,d);
        assert maps[1][2] eq W!0;
    end if;

    for n in [0..count-1] do
        stats`PointsVisited +:= 1;
        if m gt 1 and seen[n+1] then continue; end if;
        // For C!=0 the fiber is literally V12. In the overwhelmingly
        // common trivial-stabilizer case, n already encodes B.
        needVector := m gt 1 or C eq 0 or not automatic;
        if needVector then v := EF3_VectorFromKey(W,n); end if;
        if m gt 1 then
            for T in maps do
                w := v*T[1]+T[2];
                seen[EF3_Key(w)+1] := true;
            end for;
        end if;
        stats`AllOrbits +:= 1;
        if C ne 0 then
            bkey := n;
        else
            z := v*lift;
            bkey := EF3_Key(z);
        end if;
        // For canonical prefixes, the only remaining singular point is
        // C!=0,A=0,B=0. Other C=0 prefixes have A!=0 and are all smooth.
        if C ne 0 and A eq 0 and bkey eq 0 then continue; end if;
        stats`SmoothOrbits +:= 1;
        if automatic then
            ok := true;
        else
            if C ne 0 then z := v; end if;
            B := EF_FormFromVector(R,12,z);
            ok,bad := IsEverywhereMinimalF3([R|C,A,B] : Context:=Context);
            if not ok then
                assert #bad ge 1;
                if #bad eq 1 then stats`NonMinimalOne +:= 1;
                else stats`NonMinimalSeveral +:= 1; end if;
            end if;
        end if;
        if not ok then continue; end if;
        stats`K3Orbits +:= 1;
        if WriteRepresentatives then
            key := prefixPart+bkey;
            if BinaryOutput then
                for j in [0..7] do Append(~bytes,(key div 256^j) mod 256); end for;
            else
                Append(~lines,Sprintf("%o\n",key));
            end if;
            buffered +:= 1;
            if buffered eq 4096 then
                if BinaryOutput then WriteBytes(F,bytes); bytes := [Integers()|];
                else Write(F,&cat lines); lines := [Strings()|]; end if;
                buffered := 0;
            end if;
        end if;
    end for;
    if WriteRepresentatives then
        if buffered gt 0 then
            if BinaryOutput then WriteBytes(F,bytes);
            else Write(F,&cat lines); end if;
        end if;
        Flush(F);
        delete F;
    end if;
    assert stats`SmoothOrbits eq stats`K3Orbits+stats`NonMinimalOne
                                  +stats`NonMinimalSeveral;
    if complete and VerifyBurnside then
        assert stats`AllOrbits eq EF3_AffineBurnsideCount(maps);
        singular := (C ne 0 and A eq 0) select 1 else 0;
        assert stats`SmoothOrbits eq stats`AllOrbits-singular;
    end if;
    stats`CPUSeconds := Cputime(start);
    return stats;
end function;

function FibrationBenchmarkBucketsF3(prefixes : Context:=K3F3_default_context)
    // Return dictionary: <C_is_zero,d,mapCount,automatic> -> prefix indices.
    // Skip the all-singular C=A=0 stratum when forecasting the census.
    buckets := AssociativeArray();
    for i in [1..#prefixes] do
        P := prefixes[i];
        if P[1] eq 0 and P[2] eq 0 then continue; end if;
        Q,maps,S := EF3_FiberData(P);
        automatic := K3F3_PrefixAutomaticallyMinimal(P[1],P[2] : Context:=Context);
        key := <P[1] eq 0,Dimension(Q),#maps,automatic>;
        found,indices := IsDefined(buckets,key);
        if not found then indices := [Integers()|]; end if;
        Append(~indices,i);
        buckets[key] := indices;
    end for;
    return buckets;
end function;

function EF3_PrefixStem(outputStem, i)
    // Magma's Sprintf pads "%04o" with spaces, not zeros.
    digits := IntegerToString(i);
    return outputStem cat "_" cat "0"^Max(0,4-#digits) cat digits;
end function;

function CensusF3(R, indices, outputStem : WriteRepresentatives:=true,
                                         BinaryOutput:=true)
    // Each process must use disjoint indices and a distinct outputStem.
    // Existing output files cause an error. A .stats file and CSV
    // row are written only after the complete prefix returns successfully.
    prefixes := FibrationOrbitPrefixesF3(R);
    if #indices eq 0 then indices := [1..#prefixes]; end if;
    assert #indices eq #SequenceToSet(indices);
    for i in indices do assert i ge 1 and i le #prefixes; end for;
    context := K3F3_Context();
    // Check all paths before creating any outputs. This prevents stale
    // completion sidecars from surviving a restarted/truncated key file.
    EF3_RequireNewFile(outputStem cat ".csv");
    for i in indices do
        stem := EF3_PrefixStem(outputStem,i);
        EF3_RequireNewFile(stem cat ".stats");
        if WriteRepresentatives then
            EF3_RequireNewFile(stem cat (BinaryOutput select ".bin" else ".keys"));
        end if;
    end for;
    summary := Open(outputStem cat ".csv","w");
    Write(summary,"prefix,Ckey,Akey,d,maps,automatic,all,smooth,K3,nonminimal1,nonminimalSeveral,CPU_seconds,complete\n");
    result := [* *];
    for i in indices do
        P := prefixes[i];
        stem := EF3_PrefixStem(outputStem,i);
        filename := stem cat (BinaryOutput select ".bin" else ".keys");
        stats := FibrationStreamPrefixF3(P,filename :
            WriteRepresentatives:=WriteRepresentatives,BinaryOutput:=BinaryOutput,
            Context:=context);
        assert stats`Complete;
        row := Sprintf("%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
            i,EF3_Key(EF_Coefficients(P[1],4)),EF3_Key(EF_Coefficients(P[2],8)),
            stats`Dimension,stats`MapCount,stats`AutomaticMinimal,stats`AllOrbits,
            stats`SmoothOrbits,stats`K3Orbits,stats`NonMinimalOne,
            stats`NonMinimalSeveral,stats`CPUSeconds,stats`Complete);
        Write(summary,row); Flush(summary);
        printf "%o",row;
        Write(stem cat ".stats",Sprint(stats,"Magma") : Overwrite:=true);
        Append(~result,<i,stats>);
    end for;
    Flush(summary); delete summary;
    if #indices eq #prefixes then
        assert &+[ z[2]`AllOrbits : z in result ] eq 1307814237;
        assert &+[ z[2]`SmoothOrbits : z in result ] eq 1307813899;
    end if;
    return result;
end function;
