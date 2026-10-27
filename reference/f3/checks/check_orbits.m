// Run from the repository root: make reference-check/check_orbits.m
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";

procedure EF3_RunOrbitChecks()
R<s,t> := PolynomialRing(GF(3),2);
k := BaseRing(R);
clock := Cputime();
assert #EF3_BaseMatrices(k) eq 24;
for n in [0,1,2,3^13-1,3^22,3^27-1] do
    assert EF3_ModelKey(FibrationFromKeyF3(R,n)) eq n;
end for;
M := Matrix(k,2,2,[0,-1,1,0]);
N := Matrix(k,2,2,[1,1,0,1]);
g := <M,s^4+s*t^3>;
h := <N,s^3*t+t^4>;
W := [ R | s^4+s^3*t+t^4,s^8+s*t^7,s^12+t^12 ];
gh := EF3_Compose(g,h);
assert EF3_Transform(EF3_Transform(W,g[1],g[2]),h[1],h[2])
       eq EF3_Transform(W,gh[1],gh[2]);
assert FibrationDiscriminantF3(EF3_Transform(W,g[1],g[2]))
       eq EF_BaseChange(FibrationDiscriminantF3(W),g[1]);

P := FibrationOrbitPrefixesF3(R);
assert #P eq 1860;
assert #[i : i in [1..#P] | P[i][1] ne 0] eq 933;
assert #[i : i in [1..#P] | P[i][1] eq 0] eq 927;
total := FibrationStagedBurnsideCountF3(R : Prefixes := P);
assert total eq 1307814237;
singular,breakdown := SingularOrbitCountF3(R);
assert singular eq 338 and breakdown eq <23,315>;
hist,points,smoothPoints := FibrationPrefixStatisticsF3(P);
assert points eq 1493769114;
assert smoothPoints eq 1493762553;

// Check a genuinely AFFINE fiber map against the polynomial transform,
// independently of the explicit formula used to assemble FiberData.
affineChecked := false;
for prefix in P do
    if prefix[1] eq 0 then continue; end if;
    if not (&or[h[2] ne 0 : h in prefix[3]]) then continue; end if;
    Q,maps,S := EF3_FiberData(prefix);
    d := Dimension(Q);
    V12 := VectorSpace(k,13);
    for h in prefix[3] do
        origin := EF3_Transform([R|prefix[1],prefix[2],0],h[1],h[2]);
        b := Vector(k,Coordinates(Q,ReduceVector(S,
                                  V12!EF_Coefficients(origin[3],12))));
        rows := [];
        for z in Basis(Q) do
            model := [R|prefix[1],prefix[2],EF_FormFromVector(R,12,z)];
            image := EF3_Transform(model,h[1],h[2]);
            difference := V12!EF_Coefficients(image[3]-origin[3],12);
            Append(~rows,Coordinates(Q,ReduceVector(S,difference)));
        end for;
        L := Matrix(k,d,d,&cat rows);
        assert &or[T[1] eq L and T[2] eq b : T in maps];
    end for;
    affineChecked := true;
    break;
end for;
assert affineChecked;

// Enumerate a few SMALL fibers as a separate check of the scanning code.
// Both the rank-drop and ordinary j=0 branches are exercised.
tested := { Integers() | };
for prefix in P do
    if prefix[1] ne 0 then continue; end if;
    Q,maps,S := EF3_FiberData(prefix);
    d := Dimension(Q);
    if d in tested then continue; end if;
    keys,info := FibrationFiberOrbitKeysF3(prefix : GenericSmoothOnly := false);
    assert #keys eq EF3_FiberBurnsideCount(prefix);
    assert info[3] eq #keys and info[4] eq #keys;
    for key in keys do
        model := FibrationFromKeyF3(R,key);
        assert model[1] eq prefix[1] and model[2] eq prefix[2];
        assert EF3_ModelKey(model) eq key;
    end for;
    Include(~tested,d);
end for;
assert tested eq {8,9};

printf "Prefixes: 933 nonzero C + 927 zero C = %o\n",#P;
printf "All-model staged Burnside: %o\n",total;
printf "Singular orbits: %o = %o + %o\n",singular,breakdown[1],breakdown[2];
printf "Generically smooth orbits: %o\n",total-singular;
printf "All fiber scan points: %o\n",points;
printf "Scan points after skipping C=A=0: %o\n",smoothPoints;
print "Histogram <C_is_zero,dimension,distinct_maps> -> prefixes";
for key in Keys(hist) do print key,hist[key]; end for;
printf "Orbit checks passed. CPU seconds: %o\n",Cputime(clock);
end procedure;

EF3_RunOrbitChecks();
if assigned resultpath then
    Write(resultpath,"REFERENCE_CHECK_OK");
end if;
quit 0;
