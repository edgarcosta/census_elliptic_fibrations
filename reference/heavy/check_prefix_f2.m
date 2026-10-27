SetQuitOnError(true);
AttachSpec("spec");
load "reference/support/long_models.m";
load "reference/support/heavy.m";
// Compares selected full F2 fibres with permutation-group orbits and long models.
procedure test_prefix_f2()
    k := GF(2); R<s,t> := PolynomialRing(k,2);
    pm := math_base_matrices(k); V4 := math_forms(R,4); V2 := math_forms(R,2); V6 := math_forms(R,6);
    places := math_places(R,2);
    function act(W,M,q,v)
        Wb := [ math_base_change(a,M) : a in W ];
        return math_weierstrass_change(Wb,1,q*Wb[1]+q^2,q,v);
    end function;
    function sq(f) return &and[ &and[ e mod 2 eq 0 : e in Exponents(m) ] : m in Monomials(f) ]; end function;
    pre := FibrationPrefixes(R);
    targets := [ <R!0,s^6>, <R!0,s^3*t^3>, <s^2,R!0>, <s*t,R!0>, <s^2+s*t+t^2,R!0>, <s*t,s^6>, <R!0,s^5*t>, <s^2,t^6>, <s*t,s^3*t^3>, <s^2+s*t+t^2,s^6+t^6> ];
    N := 2^22; Vf := VectorSpace(k,22);
    function idx(v) return 1+&+[ Integers()!v[i]*2^(22-i) : i in [1..22] ]; end function;
    function vec(i) i -:= 1; return Vf![ (i div 2^(22-j)) mod 2 : j in [1..22] ]; end function;
    function tov(a4,a6) return Vf!(math_coefficients(a4,8) cat math_coefficients(a6,12)); end function;
    function fromv(v) return math_form(R,[v[i] : i in [1..9]],8), math_form(R,[v[9+i] : i in [1..13]],12); end function;
    Pall := Matrix(k,N,22,[ vec(i) : i in [1..N] ]);
    for tg in targets[1..heavy_target_limit(#targets)] do
        a1,a3 := Explode(tg);
        W0 := [R|a1,0,a3,0,0];
        found := 0;
        for i in [1..#pre] do
            for M in pm, q in V2 do
                W := act(W0,M,q,0);
                if W[1] eq pre[i][1] and W[3] eq pre[i][2] then found := i; break; end if;
            end for;
            if found ne 0 then break; end if;
        end for;
        assert found ne 0;
        stab := [ <M,q> : M in pm, q in V2 | W[1] eq a1 and W[3] eq a3 where W := act(W0,M,q,0) ];
        gens := [];
        for g in stab do
            M,q := Explode(g);
            z := tov(W[4],W[5]) where W := act(W0,M,q,0);
            L := Matrix(k,22,22,[ tov(W[4],W[5])-z where W := act([R|a1,0,a3,b4,b6],M,q,0) where b4,b6 := fromv(Vf.i) : i in [1..22] ]);
            P := Pall*L;
            Append(~gens, Sym(N)![ idx(P[i]+z) : i in [1..N] ]);
        end for;
        for v in math_monomials(R,6) do
            z := tov(W[4],W[5]) where W := act(W0,IdentityMatrix(k,2),R!0,v);
            Append(~gens, Sym(N)![ idx(vec(i)+z) : i in [1..N] ]);
        end for;
        PG := sub< Sym(N) | gens >;
        orbs := Orbits(PG);
        cells := AssociativeArray();
        for kd in ["smooth","regular cusp","nonregular"] do for mn in ["everywhere minimal","nonminimal somewhere"] do cells[<kd,mn>] := 0; end for; end for;
        for O in orbs do
            a4,a6 := fromv(vec(O[1]));
            W := [R|a1,0,a3,a4,a6];
            D := math_discriminant(W);
            if D ne 0 then kd := "smooth";
            elif a1 eq 0 and a3 eq 0 then kd := (sq(a4) and sq(a6)) select "nonregular" else "regular cusp";
            else kd := "nonregular"; end if;
            bad := math_bad_places(W,places,V4,V2,V6);
            cells[<kd,#bad eq 0 select "everywhere minimal" else "nonminimal somewhere">] +:= 1;
        end for;
        n,pc,rows := NumberOfFibrations(R : PrefixIndices := [found]);
        for cell in Keys(cells) do
            error if pc[cell] ne cells[cell], Sprintf(
                "prefix %o (public %o), cell %o: permutation orbits %o package %o",
                tg,found,cell,cells[cell],pc[cell]);
        end for;
        assert n eq #orbs;
    end for;
end procedure;

try
    test_prefix_f2();
catch e
    printf "MISMATCH: prefix_f2 %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
