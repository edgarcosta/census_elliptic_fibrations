SetQuitOnError(true);
AttachSpec("spec");
load "reference/support/long_models.m";
load "reference/support/heavy.m";
// Checks full long-model Burnside totals and permutation-group singular cells.
procedure test_f2_global()
    k := GF(2); R<s,t> := PolynomialRing(k,2);
    pm := math_base_matrices(k); assert #pm eq 6;
    V4 := math_forms(R,4); V2 := math_forms(R,2); V6 := math_forms(R,6); V8 := math_forms(R,8); V12 := math_forms(R,12);
    degs := [2,4,6,8,12]; n := 37; V := VectorSpace(k,n);
    function tovec(A) return V!(&cat[ math_coefficients(A[i],degs[i]) : i in [1..5] ]); end function;
    function fromvec(v) off := 0; A := [R|]; for d in degs do Append(~A,math_form(R,[v[off+i] : i in [1..d+1]],d)); off +:= d+1; end for; return A; end function;
    basis := [ fromvec(V.i) : i in [1..n] ];
    basemats := [];
    for M in pm do
        Append(~basemats, Matrix(k,n,n,[ tovec([math_base_change(a,M) : a in b]) : b in basis ]));
    end for;
    total := 0; totpre := 0; cnt := 0;
    zero := [R|0,0,0,0,0];
    V15 := VectorSpace(k,15);
    for r in V4 do for sS in V2 do for tT in V6 do
        z := tovec(math_weierstrass_change(zero,1,r,sS,tT));
        T := Matrix(k,n,n,[ tovec(math_weierstrass_change(b,1,r,sS,tT))-z : b in basis ]);
        for Bm in basemats do
            cnt +:= 1;
            L := Bm*T;
            ok := IsConsistent(L-1,-z);
            if ok then total +:= 2^(n-Rank(L-1)); end if;
            Lp := Submatrix(L,1,1,15,15); zp := V15![z[i] : i in [1..15]];
            if IsConsistent(Lp-1,-zp) then totpre +:= 2^(15-Rank(Lp-1)); end if;
        end for;
    end for; end for; end for;
    n,package_cells,rows := NumberOfFibrations(R);
    assert total eq cnt*n and totpre eq cnt*#FibrationPrefixes(R);
    places := math_places(R,2); assert #places eq 4;
    for i in [1..300] do
        a1 := Random(V2); if a1 eq 0 then continue; end if;
        A := [R|a1,Random(V4),0,0,0];
        B := math_weierstrass_change(A,1,Random(V4),Random(V2),Random(V6));
        assert math_discriminant(B) eq 0;
        ok,x0 := IsDivisibleBy(B[3],B[1]); assert ok;
        ok,y0 := IsDivisibleBy(x0^2+B[4],B[1]); assert ok;
        C := math_weierstrass_change(B,1,x0,0,y0);
        assert C[3] eq 0 and C[4] eq 0 and C[5] eq 0;
    end for;
    seen := {}; nodalreps := [];
    for a1 in V2 do if a1 eq 0 then continue; end if; for a2 in V4 do
        if <a1,a2> in seen then continue; end if;
        orb := { <math_base_change(a1,M), math_base_change(a2,M)+q*math_base_change(a1,M)+q^2> : M in pm, q in V2 };
        seen join:= orb; Append(~nodalreps,<a1,a2>);
    end for; end for;
    nodalcells := [0,0];
    for P in nodalreps do
        A := [R|P[1],P[2],0,0,0];
        bad := math_bad_places(A,places,V4,V2,V6);
        if #bad eq 0 then nodalcells[1] +:= 1; else nodalcells[2] +:= 1; end if;
        kind,mn := FibrationProperties(A); okm,cert := IsEverywhereMinimal(A);
        assert kind eq "nonregular";
        assert (mn eq "everywhere minimal") eq (#bad eq 0);
        assert cert eq <"finite",bad>;
    end for;
    for i in [1..100] do
        A := [R|0,0,0,Random(V8),Random(V12)]; q := Random(V2); v := Random(V6);
        B := math_weierstrass_change(A,1,q^2,q,v);
        assert B eq [R|0,0,0,A[4]+q^4,A[5]+q^2*A[4]+q^6+v^2];
    end for;
    N := 2^22; Vf := VectorSpace(k,22);
    function idx(v) return 1+&+[ Integers()!v[i]*2^(22-i) : i in [1..22] ]; end function;
    function vec(i) i -:= 1; return Vf![ (i div 2^(22-j)) mod 2 : j in [1..22] ]; end function;
    function tov(a4,a6) return Vf!(math_coefficients(a4,8) cat math_coefficients(a6,12)); end function;
    function fromv(v) return math_form(R,[v[i] : i in [1..9]],8), math_form(R,[v[9+i] : i in [1..13]],12); end function;
    gens := [];
    Pall := Matrix(k,N,22,[ vec(i) : i in [1..N] ]);
    Mgens := [ Matrix(k,2,2,[1,1,0,1]), Matrix(k,2,2,[0,1,1,0]) ];
    for M in Mgens do
        L := Matrix(k,22,22,[ tov(math_base_change(a4,M),math_base_change(a6,M)) where a4,a6 := fromv(Vf.i) : i in [1..22] ]);
        P := Pall*L; Append(~gens, Sym(N)![ idx(P[i]) : i in [1..N] ]);
    end for;
    for q in math_monomials(R,2) do
        z := tov(q^4,q^6);
        L := Matrix(k,22,22,[ tov(a4,a6+q^2*a4) where a4,a6 := fromv(Vf.i) : i in [1..22] ]);
        P := Pall*L; Append(~gens, Sym(N)![ idx(P[i]+z) : i in [1..N] ]);
    end for;
    for v in math_monomials(R,6) do
        z := tov(R!0,v^2);
        Append(~gens, Sym(N)![ idx(vec(i)+z) : i in [1..N] ]);
    end for;
    G := sub< Sym(N) | gens >;
    orbs := Orbits(G);
    cuspcells := [0,0,0,0];
    for O in orbs do
        a4,a6 := fromv(vec(O[1]));
        A := [R|0,0,0,a4,a6];
        sq := func< f | &and[ &and[ e mod 2 eq 0 : e in Exponents(m) ] : m in Monomials(f) ] >;
        nonreg := sq(a4) and sq(a6);
        fourth := &and[ &and[ e mod 4 eq 0 : e in Exponents(m) ] : m in Monomials(a4) ];
        bad := math_bad_places(A,places,V4,V2,V6);
        kind,mn := FibrationProperties(A); okm,cert := IsEverywhereMinimal(A);
        if nonreg then
            assert kind eq "nonregular";
            if fourth and sq(a6) then assert #bad eq 4 and cert eq <"all">; cuspcells[4] +:= 1;
            else assert cert eq <"finite",bad>; if #bad eq 0 then cuspcells[3] +:= 1; else cuspcells[4] +:= 1; end if; end if;
        else
            assert kind eq "regular cusp";
            assert &and[ TotalDegree(p) eq 1 : p in bad ];
            assert cert eq <"finite",bad>;
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
    test_f2_global();
catch e
    printf "MISMATCH: f2_global %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
