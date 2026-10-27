SetQuitOnError(true);
AttachSpec("spec");
load "reference/support/long_models.m";
load "reference/support/heavy.m";
// Checks completed-square Burnside totals and permutation-group singular cells.
procedure test_f3_global()
    k := GF(3); R<s,t> := PolynomialRing(k,2);
    pm := math_base_matrices(k); assert #pm eq 24;
    V4 := math_forms(R,4); V2 := math_forms(R,2); V6 := math_forms(R,6);
    degs3 := [4,8,12]; n3 := 27; V3 := VectorSpace(k,n3);
    function tov3(W) return V3!(&cat[ math_coefficients(W[i],degs3[i]) : i in [1..3] ]); end function;
    function fromv3(v) off := 0; W := [R|]; for d in degs3 do Append(~W,math_form(R,[v[off+i] : i in [1..d+1]],d)); off +:= d+1; end for; return W; end function;
    function act3(W,M,u,q)
        C := u^2*math_base_change(W[1],M); A := u^4*math_base_change(W[2],M); B := u^6*math_base_change(W[3],M);
        return [R| C, A-C*q, B+A*q+C*q^2+q^3 ];
    end function;
    for i in [1..50] do
        W := [ Random(V4), math_form(R,Random(VectorSpace(k,9)),8), math_form(R,Random(VectorSpace(k,13)),12) ]; M := Random(pm); q := Random(V4); u := Random([k|1,2]);
        Lb := [ math_base_change(a,M) : a in [R|0,W[1],0,W[2],W[3]] ];
        L2 := math_weierstrass_change(Lb,u,q,0,0);
        assert [L2[2],L2[4],L2[5]] eq act3(W,M,u,q) and L2[1] eq 0 and L2[3] eq 0;
    end for;
    b3 := [ fromv3(V3.i) : i in [1..n3] ];
    total := 0; totalpre := 0; cnt := 0;
    V14 := VectorSpace(k,14);
    for M in pm do for u in [k|1,2] do for q in V4 do
        cnt +:= 1;
        z := tov3(act3([R|0,0,0],M,u,q));
        L := Matrix(k,n3,n3,[ tov3(act3(b,M,u,q))-z : b in b3 ]);
        ok,x := IsConsistent(L-1, -z);
        if ok then total +:= 3^(n3-Rank(L-1)); end if;
        Lp := Submatrix(L,1,1,14,14); zp := V14![z[i] : i in [1..14]];
        okp := IsConsistent(Lp-1,-zp);
        if okp then totalpre +:= 3^(14-Rank(Lp-1)); end if;
    end for; end for; end for;
    n,package_cells,rows := NumberOfFibrations(R);
    assert total eq cnt*n and totalpre eq cnt*#FibrationPrefixes(R);
    assert cnt eq 24*2*243;
    for i in [1..200] do
        sS := Random(V2); tT := Random(V6); r := Random(V4);
        W := [R|0,Random(V4),0,math_form(R,Random(VectorSpace(k,9)),8),math_form(R,Random(VectorSpace(k,13)),12)];
        L := math_weierstrass_change(W,1,r,sS,tT);
        assert (L[1] eq 0 and L[3] eq 0) eq (sS eq 0 and tT eq 0);
    end for;
    found := 0; tested := 0;
    for i in [1..3000] do
        C := Random(V4); if C eq 0 then continue; end if;
        A := math_form(R,Random(VectorSpace(k,9)),8);
        ok,B := IsDivisibleBy(A^2*(C^2-A),C^3);
        if not ok then continue; end if;
        tested +:= 1;
        W := [R|C,A,B]; assert C^2*A^2-A^3-C^3*B eq 0;
        okq,q := IsDivisibleBy(A,C);
        assert okq and TotalDegree(q) le 4;
        W2 := act3(W,IdentityMatrix(k,2),1,q);
        assert W2[2] eq 0 and W2[3] eq 0 and W2[1] eq C;
        found +:= 1;
    end for;
    quartics := [ C : C in V4 | C ne 0 ];
    seen := {}; nodalreps := [R|];
    for C in quartics do
        if C in seen then continue; end if;
        orb := { math_base_change(C,M) : M in pm }; seen join:= orb; Append(~nodalreps,C);
    end for;
    places := math_places(R,2); assert #places eq 7;
    nodalcells := [0,0];
    for C in nodalreps do
        A := [R|0,C,0,0,0];
        bad := math_bad_places(A,places,V4,V2,V6);
        if #bad eq 0 then nodalcells[1] +:= 1; else nodalcells[2] +:= 1; end if;
        kind,mn := FibrationProperties(A); okm,cert := IsEverywhereMinimal(A);
        assert kind eq "nonregular";
        assert (mn eq "everywhere minimal") eq (#bad eq 0);
        assert cert eq <"finite",bad>;
    end for;
    V12 := VectorSpace(k,13);
    N := 3^13;
    function idx(v) return 1+&+[ Integers()!v[i]*3^(13-i) : i in [1..13] ]; end function;
    function vec(i) i -:= 1; return V12![ (i div 3^(13-j)) mod 3 : j in [1..13] ]; end function;
    gens := [];
    B12 := math_monomials(R,12);
    Mgens := [ Matrix(k,2,2,[1,1,0,1]), Matrix(k,2,2,[0,1,1,0]), Matrix(k,2,2,[1,0,0,2]) ];
    for M in Mgens do
        L := Matrix(k,13,13,[ math_coefficients(math_base_change(B12[i],M),12) : i in [1..13] ]);
        P := Matrix(k,N,13,[ vec(i) : i in [1..N] ])*L;
        Append(~gens, Sym(N)![ idx(P[i]) : i in [1..N] ]);
    end for;
    for m in math_monomials(R,4) do
        c := V12!math_coefficients(m^3,12);
        Append(~gens, Sym(N)![ idx(vec(i)+c) : i in [1..N] ]);
    end for;
    G := sub< Sym(N) | gens >;
    orbs := Orbits(G);
    cuspcells := [0,0,0,0];
    for O in orbs do
        B := math_form(R,vec(O[1]),12);
        A := [R|0,0,0,0,B];
        iscube := &and[ &and[ e mod 3 eq 0 : e in Exponents(m) ] : m in Monomials(B) ];
        bad := math_bad_places(A,places,V4,V2,V6);
        kind,mn := FibrationProperties(A); okm,cert := IsEverywhereMinimal(A);
        if iscube then
            assert #bad eq 7;
            assert kind eq "nonregular" and cert eq <"all">;
            cuspcells[4] +:= 1;
        else
            assert kind eq "regular cusp";
            assert &and[ TotalDegree(p) eq 1 : p in bad ];
            assert cert eq <"finite",bad>;
            assert (mn eq "everywhere minimal") eq (#bad eq 0);
            if #bad eq 0 then cuspcells[1] +:= 1; else cuspcells[2] +:= 1; end if;
        end if;
    end for;
    assert cuspcells[1] eq package_cells[<"regular cusp","everywhere minimal">];
    assert cuspcells[2] eq package_cells[<"regular cusp","nonminimal somewhere">];
    assert cuspcells[3]+nodalcells[1] eq package_cells[<"nonregular","everywhere minimal">];
    assert cuspcells[4]+nodalcells[2] eq package_cells[<"nonregular","nonminimal somewhere">];
    assert #orbs eq #FibrationRepresentativeKeys(R,Index(FibrationPrefixes(R),<R!0,R!0>));
end procedure;

try
    test_f3_global();
catch e
    printf "MISMATCH: f3_global %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
