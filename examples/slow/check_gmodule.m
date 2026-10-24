// Compares ef3_gmodule.m with the scan-based code in ef3_orbits.m,
// k3_minimality_f3.m and ef3_stream.m.
// Run from the repository root: magma -b filename:=examples/slow/check_gmodule.m exitsignal:='' run_tests.m

procedure EF3G_RunChecks()
clock := Cputime();
k := GF(3);
R<s,t> := PolynomialRing(k,2);
G := GL(2,k);
for d in [4,8,12] do
    rho := EF_FormAction(R,G,d);
    for g in G do
        assert rho(g) eq EF_FormRep(R,d,Matrix(g));
        for h in G do assert rho(g*h) eq rho(g)*rho(h); end for;
    end for;
end for;

old := FibrationOrbitPrefixesF3(R);
new := FibrationOrbitPrefixesF3G(R);
assert #new eq 1860 and #[ P : P in new | P[1] ne 0 ] eq 933;

// A pair (C,A) up to the group, for C a least-key orbit representative:
// the least key of A modulo C*V4 over the stabilizer of C.
V8 := VectorSpace(k,9);
B4 := EF_BinaryBasis(R,4);
base := EF3_BaseMatrices(k);
function Canonical(C, A)
    S := sub< V8 | [ V8!EF_Coefficients(C*m,8) : m in B4 ] >;
    keys := [ EF3_Key(ReduceVector(S,V8!EF_Coefficients(EF_BaseChange(A,M),8)))
              : M in base | EF_BaseChange(C,M) eq C ];
    return <EF3_Key(EF_Coefficients(C,4)),Min(keys)>;
end function;
oldIndex := AssociativeArray();
for i in [1..#old] do oldIndex[Canonical(old[i][1],old[i][2])] := i; end for;
match := [];
for P in new do
    i := oldIndex[Canonical(P[1],P[2])];
    assert 2*#old[i][3] eq #P[3];
    Append(~match, i);
end for;
assert #SequenceToSet(match) eq #old;
print "Prefixes match the scan-based ones, with the same stabilizer orders.";

fibres := [ FibrationFibreF3G(P) : P in new ];
counts := [ EF_FibreCounts(F) : F in fibres ];
totals := < &+[ c[1] : c in counts ], &+[ c[2] : c in counts ], &+[ c[3] : c in counts ] >;
assert totals[1] eq 1307814237 and totals[2] eq 1307813899;
for j in [1..#new] do
    assert counts[j][1] eq EF3_FiberBurnsideCount(old[match[j]]);
end for;
printf "Fibre counts match prefix by prefix; totals %o.\n", totals;

// Elliptic K3 counts against enumeration with the local test: every pair
// with C = 0, and for C != 0 one pair that needs the local test for each
// nontrivial stabilizer order on the fibre.
compared := 0;
orders := {};
for j in [1..#new] do
    P := new[j];
    if P[1] eq 0 and P[2] eq 0 then continue; end if;
    F := fibres[j];
    d := F`Dimension;
    E := F`Group;
    automatic := K3F3_PrefixAutomaticallyMinimal(P[1],P[2]);
    if P[1] ne 0 and (automatic or #E eq 1 or #E in orders) then continue; end if;
    if P[1] ne 0 then Include(~orders, #E); end if;
    st := FibrationStreamPrefixF3(old[match[j]],"unused" :
                                  WriteRepresentatives := false);
    assert st`K3Orbits eq counts[j][3];
    compared +:= 1;
end for;
printf "Elliptic K3 counts match enumeration on %o pairs.\n", compared;

// Representatives and the non-minimal set against the local test.
kinds := {};
for j in [1..#new] do
    P := new[j];
    if P[1] eq 0 and P[2] eq 0 then continue; end if;
    if K3F3_PrefixAutomaticallyMinimal(P[1],P[2]) then continue; end if;
    F := fibres[j];
    d := F`Dimension;
    E := F`Group;
    kind := <P[1] eq 0, d, #E gt 1>;
    if kind in kinds then continue; end if;
    Include(~kinds, kind);
    reps := EF_FibreRepresentatives(F);
    assert #reps eq counts[j][3];
    // With a trivial stabilizer every point is its own orbit; sample those.
    if #E eq 1 then reps := [ Random(reps) : i in [1..20000] ]; end if;
    for key in reps do
        assert IsEverywhereMinimalF3(FibrationFromKeyF3(R,key));
    end for;
    W := VectorSpace(k,d+1);
    canon := { Min([ EF3_Key(Eltseq(W!(Eltseq(F`FromKey(key)) cat [1])*g))
                     : g in E ]) : key in reps };
    assert #E eq 1 or #canon eq #reps;
    bad := Set(EF_NonMinimalPoints(F));
    for v in bad join { Random(VectorSpace(k,d)) : i in [1..500] } do
        excluded := v in bad or EF_InSingular(F,v);
        assert excluded eq not IsEverywhereMinimalF3(FibrationFromKeyF3(R,F`Key(v)));
    end for;
    printf "Pair %o, kind %o: %o representatives checked.\n", j, kind, #reps;
end for;
// Representatives over a fibre whose stabilizer fixes no point.
for j in [1..#new] do
    P := new[j];
    F := fibres[j];
    d := F`Dimension;
    E := F`Group;
    if #E eq 1 then continue; end if;
    X := [ Matrix(E.i) : i in [1..Ngens(E)] ];
    if IsConsistent(HorizontalJoin([ Submatrix(x,1,1,d,d)-IdentityMatrix(k,d)
                                     : x in X ]),
                    Vector(k, &cat[ [ -x[d+1,i] : i in [1..d] ] : x in X ]))
    then continue; end if;
    reps := EF_FibreRepresentatives(F);
    assert #reps eq counts[j][3];
    W := VectorSpace(k,d+1);
    canon := { Min([ EF3_Key(Eltseq(W!(Eltseq(F`FromKey(key)) cat [1])*g))
                     : g in E ]) : key in reps };
    assert #canon eq #reps;
    printf "Pair %o, no fixed point: %o representatives checked.\n", j, #reps;
    break;
end for;
printf "All checks passed. CPU seconds: %o\n", Cputime(clock);
end procedure;

EF3G_RunChecks();
