AttachSpec("spec");
load "reference/f2/elliptic_fibration_orbits_f2.m";
load "reference/f2/k3_minimality_f2.m";

// Checks all F2 counts and smooth minimal representatives by reference orbits.
procedure test_f2_scan_comparison()
    k := GF(2);
    R<s,t> := PolynomialRing(k,2);
    old := FibrationOrbitPrefixesF2(R);
    new := FibrationPrefixes(R);
    assert #old eq 88 and #new eq 88;
    pairs := [ <Matrix(g),EF2FormFromMask(R,2,n)> : g in GL(2,k), n in [0..7] ];
    function canonical_prefix(a1,a3)
        return Min([ EF2BitKey(EF2Coefficients(B[1],2) cat EF2Coefficients(B[3],6))
            where B is TransformNormalizedF2([R | a1,0,a3,0,0],h[1],h[2],R!0)
            : h in pairs ]);
    end function;
    index := AssociativeArray();
    for j in [1..#new] do index[canonical_prefix(new[j][1],new[j][2])] := j; end for;
    assert #Keys(index) eq #old;
    n, cells, rows := NumberOfFibrations(R);
    totals := [0,0,0];
    for P in old do
        j := index[canonical_prefix(P[1],P[2])];
        Q := new[j];
        transport := [ h : h in pairs | B[1] eq P[1] and B[3] eq P[2]
            where B is TransformNormalizedF2([R | Q[1],0,Q[2],0,0],h[1],h[2],R!0) ];
        assert not IsEmpty(transport);
        h := transport[1];
        keys, info := FibrationFiberOrbitKeysF2(P);
        minimal := [ key : key in keys | IsEverywhereMinimalF2(FibrationFromKeyF2(R,key)) ];
        counts := rows[j][2];
        c := < &+[ counts[cell] : cell in Keys(counts) ],
            counts[<"smooth","everywhere minimal">]+counts[<"smooth","nonminimal somewhere">],
            counts[<"smooth","everywhere minimal">] >;
        assert c eq <info[3],info[4],#minimal>;
        for i in [1..3] do totals[i] +:= c[i]; end for;

        V := VectorSpace(k,22);
        S := sub< V | [ V!EF2PairVector(P[1]*r,r^2+P[2]*r)
            : r in EF2BinaryBasis(R,6) ] >;
        C := sub< V | [ ReduceVector(S,v) : v in Basis(V) ] >;
        stabilizer := [ g : g in pairs | B[1] eq P[1] and B[3] eq P[2]
            where B is TransformNormalizedF2([R | P[1],0,P[2],0,0],g[1],g[2],R!0) ];
        maps := [* *];
        for g in stabilizer do
            origin := TransformNormalizedF2([R | P[1],0,P[2],0,0],g[1],g[2],R!0);
            b := C!ReduceVector(S,V!EF2PairVector(origin[4],origin[5]));
            images := [];
            for v in Basis(C) do
                a4,a6 := EF2PairFromVector(R,v);
                W := TransformNormalizedF2([R | P[1],0,P[2],a4,a6],g[1],g[2],R!0);
                Append(~images,Coordinates(C,C!ReduceVector(S,V!EF2PairVector(W[4],W[5]))-b));
            end for;
            Append(~maps,<Matrix(k,images),Vector(k,Coordinates(C,b))>);
        end for;
        function orbit_key(W)
            v := Vector(k,Coordinates(C,C!ReduceVector(S,V!EF2PairVector(W[4],W[5]))));
            return Min([ EF2BitKey(Eltseq(v*g[1]+g[2])) : g in maps ]);
        end function;
        expected := { orbit_key(FibrationFromKeyF2(R,key)) : key in minimal };
        reps := FibrationRepresentativeKeys(R,j : GenericFibreKinds := ["smooth"],
            Minimalities := ["everywhere minimal"]);
        assert #reps eq #minimal and #expected eq #minimal;
        actual := { orbit_key(TransformNormalizedF2(aInvariantsFromKey(R,key),h[1],h[2],R!0))
            : key in reps };
        assert actual eq expected;
    end for;
    assert n eq totals[1] and totals eq [1402468,1401732,1397598];
end procedure;

test_f2_scan_comparison();
Write(resultpath,"REFERENCE_CHECK_OK");
quit 0;
