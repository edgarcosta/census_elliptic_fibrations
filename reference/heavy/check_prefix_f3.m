SetQuitOnError(true);
AttachSpec("spec");
load "reference/support/long_models.m";
load "reference/support/heavy.m";
// Compares selected full F3 fibres with permutation-group orbits and long models.
procedure test_prefix_f3()
    k := GF(3); R<s,t> := PolynomialRing(k,2);
    pm := math_base_matrices(k); V4 := math_forms(R,4); V2 := math_forms(R,2); V6 := math_forms(R,6);
    places := math_places(R,2);
    G := [ <M,q> : M in pm, q in V4 ]; assert #G eq 5832;
    function act(W,g) M,q := Explode(g); C := math_base_change(W[1],M); A := math_base_change(W[2],M); B := math_base_change(W[3],M); return [R|C,A-C*q,B+A*q+C*q^2+q^3]; end function;
    pre := FibrationPrefixes(R);
    targets := [ <s^4,R!0>, <s^2*t^2,R!0>, <(s^2+t^2)^2,R!0>, <s^3*t,R!0>, <s^4,t^8>, <s^4+t^4,R!0>, <s^2*t^2,s^8+t^8>, <s^4,s^4*t^4>,
                 <R!0,s^8>, <R!0,s^4*t^4>, <R!0,(s^2+t^2)^4>, <R!0,s^7*t>, <R!0,s^4*t^4+t^8>, <s*t*(s^2+t^2),R!0>, <s^2*t^2,s^6*t^2> ];
    N := 3^13; V12 := VectorSpace(k,13);
    function idx(v) return 1+&+[ Integers()!v[i]*3^(13-i) : i in [1..13] ]; end function;
    function vec(i) i -:= 1; return V12![ (i div 3^(13-j)) mod 3 : j in [1..13] ]; end function;
    Pall := Matrix(k,N,13,[ vec(i) : i in [1..N] ]);
    B12 := math_monomials(R,12);
    for tg in targets[1..heavy_target_limit(#targets)] do
        C,A := Explode(tg);
        images := [ act([R|C,A,0],g) : g in G ];
        orbitset := { <W[1],W[2]> : W in images };
        found := [ i : i in [1..#pre] | <pre[i][1],pre[i][2]> in orbitset ];
        assert #found eq 1; found := found[1];
        stab := [ G[j] : j in [1..#G] | images[j][1] eq C and images[j][2] eq A ];
        gens := [];
        for g in stab do
            z := V12!math_coefficients(act([R|C,A,0],g)[3],12);
            L := Matrix(k,13,13,[ V12!math_coefficients(act([R|C,A,B12[i]],g)[3],12)-z : i in [1..13] ]);
            P := Pall*L;
            Append(~gens, Sym(N)![ idx(P[i]+z) : i in [1..N] ]);
        end for;
        PG := sub< Sym(N) | gens >;
        orbs := Orbits(PG);
        cells := AssociativeArray();
        for kd in ["smooth","regular cusp","nonregular"] do for mn in ["everywhere minimal","nonminimal somewhere"] do cells[<kd,mn>] := 0; end for; end for;
        for O in orbs do
            B := math_form(R,vec(O[1]),12);
            D := C^2*A^2-A^3-C^3*B;
            if D ne 0 then kd := "smooth";
            elif C eq 0 and A eq 0 then kd := (&and[ &and[ e mod 3 eq 0 : e in Exponents(m) ] : m in Monomials(B) ]) select "nonregular" else "regular cusp";
            else kd := "nonregular"; end if;
            bad := math_bad_places([R|0,C,0,A,B],places,V4,V2,V6);
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
    test_prefix_f3();
catch e
    printf "MISMATCH: prefix_f3 %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
