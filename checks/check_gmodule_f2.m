// Compares ef2_gmodule.m with the scan-based code in
// elliptic_fibration_orbits_f2.m and k3_minimality_f2.m on all 88 prefixes:
// orbit counts, and the K3 representatives up to isomorphism.
// Run from the repository root: magma -b checks/check_gmodule_f2.m
load "elliptic_fibration_orbits_f2.m";
load "k3_minimality_f2.m";
load "ef_common.m";
load "ef_orbits.m";
load "ef2_gmodule.m";

procedure EF2G_RunChecks()
clock := Cputime();
k := GF(2);
R<s,t> := PolynomialRing(k,2);
G := GL(2,k);
for d in [2,6,8,12] do
    rho := EF_FormAction(R,G,d);
    for g, h in G do
        assert rho(g) eq EF_FormRep(R,d,Matrix(g)) and rho(g*h) eq rho(g)*rho(h);
    end for;
end for;

old := FibrationOrbitPrefixesF2(R);
new := FibrationOrbitPrefixesF2G(R);
assert #old eq 88 and #new eq 88;

// A prefix up to the group: the least key of (a1,a3) over all <M,p>, and a
// <M,p> reaching it.
pairs := [ <Matrix(g),EF2FormFromMask(R,2,n)> : g in G, n in [0..7] ];
function Canonical(a1, a3)
    best := 2^10;
    for h in pairs do
        B := TransformNormalizedF2([ R | a1,0,a3,0,0 ], h[1], h[2], R!0);
        key := EF2BitKey(EF_Coefficients(B[1],2) cat EF_Coefficients(B[3],6));
        if key lt best then best := key; arg := h; end if;
    end for;
    return best, arg;
end function;
newIndex := AssociativeArray();
for j in [1..#new] do newIndex[Canonical(new[j][1],new[j][2])] := j; end for;

totals := [0,0,0];
for P in old do
    keyP, toCanonical := Canonical(P[1],P[2]);
    j := newIndex[keyP];
    Q := new[j];
    F := FibrationFibreF2G(Q);
    c := EF_FibreCounts(F);
    keys, info := FibrationFiberOrbitKeysF2(P);
    k3 := [ key : key in keys | IsEverywhereMinimalF2(FibrationFromKeyF2(R,key)) ];
    assert c eq <info[3], info[4], #k3>;
    for i in [1..3] do totals[i] +:= c[i]; end for;

    // Carry the old prefix to the new one: to the canonical pair, then back.
    found := false;
    for h in pairs do
        B := TransformNormalizedF2(TransformNormalizedF2([ R | P[1],0,P[2],0,0 ],
                 toCanonical[1], toCanonical[2], R!0), h[1], h[2], R!0);
        if B[1] eq Q[1] and B[3] eq Q[2] then found := true; back := h; break; end if;
    end for;
    assert found;
    d := F`Dimension;
    E := F`Group;
    W := VectorSpace(k,d+1);
    function OrbitKey(v)
        return Min([ EF_PrimeFieldKey(Eltseq(W!(Eltseq(v) cat [1])*g),2) : g in E ]);
    end function;
    oldOrbits := { OrbitKey(F`FromKey(EF2BitKey(&cat[ EF_Coefficients(B[i],[2,0,6,8,12][i])
                     : i in [1,3,4,5] ])))
                   where B is TransformNormalizedF2(TransformNormalizedF2(
                       FibrationFromKeyF2(R,key), toCanonical[1], toCanonical[2], R!0),
                       back[1], back[2], R!0)
                   : key in k3 };
    reps := EF_FibreRepresentatives(F);
    assert #reps eq #k3;
    assert { OrbitKey(F`FromKey(key)) : key in reps } eq oldOrbits;
end for;
assert totals eq [1402468,1401732,1397598];
printf "Counts and K3 representatives match on all 88 prefixes; totals %o.\n", totals;
printf "All checks passed. CPU seconds: %o\n", Cputime(clock);
end procedure;

EF2G_RunChecks();
quit;
