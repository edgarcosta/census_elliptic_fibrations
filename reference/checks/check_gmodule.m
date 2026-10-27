AttachSpec("spec");
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/f3/ef3_stream.m";

// Checks public F3 counts and representative membership using reference maps.
procedure test_f3_scan_comparison()
    k := GF(3);
    R<s,t> := PolynomialRing(k,2);
    old := FibrationOrbitPrefixesF3(R);
    new := FibrationPrefixes(R);
    assert #old eq 1860 and #new eq 1860;
    V8 := VectorSpace(k,9);
    B4 := EF_BinaryBasis(R,4);
    base := EF3_BaseMatrices(k);
    function canonical_prefix(C,A)
        candidates := [];
        for M in base do
            D := EF_BaseChange(C,M);
            S := sub< V8 | [ V8!EF_Coefficients(D*r,8) : r in B4 ] >;
            a := ReduceVector(S,V8!EF_Coefficients(EF_BaseChange(A,M),8));
            Append(~candidates,<EF3_Key(EF_Coefficients(D,4)),EF3_Key(a)>);
        end for;
        return Min(candidates);
    end function;
    index := AssociativeArray();
    for i in [1..#old] do index[canonical_prefix(old[i][1],old[i][2])] := i; end for;
    assert #Keys(index) eq #old;
    match := [ index[canonical_prefix(P[1],P[2])] : P in new ];
    assert #SequenceToSet(match) eq #old;
    n, cells, rows := NumberOfFibrations(R);
    counts := [ < &+[ c[cell] : cell in Keys(c) ],
        c[<"smooth","everywhere minimal">]+c[<"smooth","nonminimal somewhere">],
        c[<"smooth","everywhere minimal">] > where c is row[2] : row in rows ];
    assert n eq 1307814237;
    assert &+[ c[2] : c in counts ] eq 1307813899;
    assert &+[ c[3] : c in counts ] eq 1307723786;
    for j in [1..#new] do
        assert counts[j][1] eq EF3_FiberBurnsideCount(old[match[j]]);
    end for;

    // Match every C=0 smooth fibre and one nonautomatic C!=0 per map order.
    orders := {};
    compared := 0;
    kinds := {};
    no_fixed_point := false;
    for j in [1..#new] do
        P := old[match[j]];
        if P[1] eq 0 and P[2] eq 0 then continue; end if;
        Q,maps,S := EF3_FiberData(P);
        d := Dimension(Q);
        automatic := K3F3_PrefixAutomaticallyMinimal(P[1],P[2]);
        if P[1] eq 0 or (not automatic and #maps gt 1 and #maps notin orders) then
            stats := FibrationStreamPrefixF3(P,"unused" : WriteRepresentatives := false);
            assert stats`AllOrbits eq counts[j][1];
            assert stats`SmoothOrbits eq counts[j][2];
            assert stats`K3Orbits eq counts[j][3];
            if P[1] ne 0 then Include(~orders,#maps); end if;
            compared +:= 1;
        end if;
        fixed := IsConsistent(HorizontalJoin([ g[1]-IdentityMatrix(k,d) : g in maps ]),
            Vector(k,&cat[ Eltseq(-g[2]) : g in maps ]));
        kind := <P[1] eq 0,d,#maps gt 1>;
        check_kind := not automatic and kind notin kinds;
        check_affine := not fixed and not no_fixed_point;
        if not check_kind and not check_affine then continue; end if;
        if check_kind then Include(~kinds,kind); end if;
        if check_affine then no_fixed_point := true; end if;

        // Carry public representatives to the matched reference prefix.
        found := false;
        for M in base do
            if EF_BaseChange(new[j][1],M) ne P[1] then continue; end if;
            A := EF_BaseChange(new[j][2],M);
            if P[1] eq 0 then
                if A ne P[2] then continue; end if;
                q := R!0;
            else
                divisible,q := IsDivisibleBy(A-P[2],P[1]);
                if not divisible then continue; end if;
            end if;
            transport := <M,q>;
            found := true;
            break;
        end for;
        assert found;
        function orbit_key(W)
            v := Vector(k,Coordinates(Q,Q!ReduceVector(S,
                Vector(k,EF_Coefficients(W[3],12)))));
            return Min([ EF3_Key(v*g[1]+g[2]) : g in maps ]);
        end function;
        reps := FibrationRepresentativeKeys(R,j : GenericFibreKinds := ["smooth"],
            Minimalities := ["everywhere minimal"]);
        assert not IsEmpty(reps) and #reps eq counts[j][3];
        // Trivial fibres can have millions of points; check a bounded sample.
        if #maps eq 1 then reps := [ reps[i] : i in [1..Min([#reps,20000])] ]; end if;
        actual := {};
        for key in reps do
            W := EF3_NormalizeLongModel(aInvariantsFromKey(R,key));
            W := EF3_Transform(W,transport[1],transport[2]);
            assert W[1] eq P[1] and W[2] eq P[2];
            assert IsEverywhereMinimalF3(W);
            Include(~actual,orbit_key(W));
        end for;
        assert #actual eq #reps;
        if #maps gt 1 then
            keys,info := FibrationFiberOrbitKeysF3(P);
            expected := { orbit_key(W) : key in keys | IsEverywhereMinimalF3(W)
                where W is FibrationFromKeyF3(R,key) };
            assert actual eq expected;
        end if;
    end for;
    assert compared eq 930 and not IsEmpty(kinds) and no_fixed_point;
end procedure;

test_f3_scan_comparison();
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
