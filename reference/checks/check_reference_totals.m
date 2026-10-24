AttachSpec("spec");
load "reference/f2/elliptic_fibration_orbits_f2.m";
load "reference/f2/k3_minimality_f2.m";
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/support/classification.m";
load "reference/support/singular.m";

function reference_f2_burnside(P)
    a1,a3 := Explode(P);
    R := Parent(a1); k := BaseRing(R);
    V := VectorSpace(k,22);
    S := sub< V | [ V!EF2PairVector(a1*r,r^2+a3*r) : r in EF2BinaryBasis(R,6) ] >;
    Q := sub< V | [ ReduceVector(S,v) : v in Basis(V) ] >;
    d := Dimension(Q);
    fixed := 0; order := 0;
    for g in GL(2,k), n in [0..7] do
        M := Matrix(g); q := EF2FormFromMask(R,2,n);
        B := TransformNormalizedF2([ R | a1,0,a3,0,0 ],M,q,R!0);
        if B[1] ne a1 or B[3] ne a3 then continue; end if;
        b := Vector(k,Coordinates(Q,Q!ReduceVector(S,V!EF2PairVector(B[4],B[5]))));
        rows := [];
        for v in Basis(Q) do
            a4,a6 := EF2PairFromVector(R,v);
            B := TransformNormalizedF2([ R | a1,0,a3,a4,a6 ],M,q,R!0);
            Append(~rows,Eltseq(Vector(k,Coordinates(Q,Q!ReduceVector(S,
                V!EF2PairVector(B[4],B[5]))))-b));
        end for;
        A := Matrix(k,rows)-IdentityMatrix(k,d);
        if IsConsistent(A,-b) then fixed +:= 2^(d-Rank(A)); end if;
        order +:= 1;
    end for;
    assert order gt 0 and fixed mod order eq 0;
    return fixed div order;
end function;

function matching_prefix(P,prefixes)
    R := Parent(P[1]); p := Characteristic(R);
    for i in [1..#prefixes] do
        Q := prefixes[i];
        for g in GL(2,GF(p)) do
            M := Matrix(g);
            if p eq 2 then
                if EF2BaseChange(Q[1],M) ne P[1] then continue; end if;
                for n in [0..7] do
                    B := TransformNormalizedF2([ R | Q[1],0,Q[2],0,0 ],M,
                        EF2FormFromMask(R,2,n),R!0);
                    if B[1] eq P[1] and B[3] eq P[2] then return i; end if;
                end for;
            else
                if EF_BaseChange(Q[1],M) ne P[1] then continue; end if;
                A := EF_BaseChange(Q[2],M)-P[2];
                if P[1] eq 0 then
                    if A eq 0 then return i; end if;
                elif IsDivisibleBy(A,P[1]) then
                    return i;
                end if;
            end if;
        end for;
    end for;
    error "No public prefix matches the reference prefix",P;
end function;

// Complete selected fibres detect a wrong split even if their totals agree.
// F2 samples each a1; F3 samples each C=0 map order and minimality shortcut.
// check_gmodule_f2.m recounts all F2 fibres; check_gmodule.m covers more F3
// fibres. The F3 smooth minimal/nonminimal split is sampled, not exhaustive.
procedure test_reference_totals()
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        models := reference_singular_models(R);
        singular := reference_singular_cells(models);
        old := p eq 2 select FibrationOrbitPrefixesF2(R) else FibrationOrbitPrefixesF3(R);
        all_count := p eq 2 select &+[ reference_f2_burnside(P) : P in old ]
            else &+[ EF3_FiberBurnsideCount(P) : P in old ];
        smooth_count := all_count-&+singular;
        assert all_count eq (p eq 2 select 1402468 else 1307814237);
        seen := {};
        sample := [0,0];
        samples := [* *];
        for P in old do
            if P[1] eq 0 and P[2] eq 0 then continue; end if;
            if p eq 2 then
                tag := P[1];
                if tag in seen then continue; end if;
                keys := FibrationFiberOrbitKeysF2(P);
                good := #[ key : key in keys | IsEverywhereMinimalF2(FibrationFromKeyF2(R,key)) ];
            else
                if P[1] ne 0 then continue; end if;
                Q,maps := EF3_FiberData(P);
                tag := <#maps,K3F3_PrefixAutomaticallyMinimal(P[1],P[2])>;
                if tag in seen then continue; end if;
                keys := FibrationFiberOrbitKeysF3(P);
                good := #[ key : key in keys | IsEverywhereMinimalF3(FibrationFromKeyF3(R,key)) ];
            end if;
            Include(~seen,tag);
            expected := [good,#keys-good];
            Append(~samples,<P,expected>);
            for i in [1..2] do sample[i] +:= expected[i]; end for;
        end for;
        assert #samples eq (p eq 2 select 4 else 13);
        assert sample[1] gt 0 and sample[2] gt 0;
        n,cells,rows := NumberOfFibrations(R);
        assert n eq all_count;
        assert cells[<"smooth","everywhere minimal">]+cells[<"smooth","nonminimal somewhere">] eq smooth_count;
        assert [ cells[<kind,minimum>] : minimum in ["everywhere minimal","nonminimal somewhere"],
            kind in ["regular cusp","nonregular"] ] eq singular;
        prefixes := FibrationPrefixes(R);
        for item in samples do
            j := matching_prefix(item[1],prefixes);
            actual := rows[j][2];
            assert actual[<"smooth","everywhere minimal">] eq item[2][1];
            assert actual[<"smooth","nonminimal somewhere">] eq item[2][2];
        end for;
        printf "F%o reference totals: all=%o smooth=%o; singular cells=%o\n",
            p,all_count,smooth_count,singular;
        printf "F%o complete-fibre smooth sample: prefixes=%o minimal=%o nonminimal=%o\n",
            p,#samples,sample[1],sample[2];
    end for;
end procedure;

test_reference_totals();
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
