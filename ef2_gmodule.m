// The census over GF(2) for ef_orbits.m: models
// y^2 + a1*x*y + a3*y = x^3 + a4*x + a6, a2 having been removed by
// x -> x+a2, with [a1,a3,a4,a6] of degrees 2, 6, 8, 12. Load ef_common.m and
// ef_orbits.m first. A prefix is <a1,a3,H> with H the pairs <M,p> in
// GL(2,2) x V2 for which M followed by x -> x+p^2+a1*p, y -> y+p*x fixes
// a1 and a3. Keys are those of elliptic_fibration_orbits_f2.m.

function FibrationOrbitPrefixesF2G(R)
    k := BaseRing(R);
    assert #k eq 2 and Rank(R) eq 2;
    G := GL(2,k);
    rho2 := EF_FormAction(R,G,2);
    rho6 := EF_FormAction(R,G,6);
    G2 := Image(rho2);
    B2 := EF_BinaryBasis(R,2);
    V2 := VectorSpace(k,3);
    prefixes := [* *];
    for c in EF_OrbitRepresentatives(G2) do
        a1 := EF_FormFromVector(R,2,c);
        H1 := Stabiliser(G2,c) @@ rho2;
        // The normalized y -> y+p*x shifts of a3 are GF(2)-linear.
        shifts := [ EF_WeierstrassChange([ R | a1,0,0,0,0 ],
                        k!1,m^2+a1*m,m,R!0)[3] : m in B2 ];
        P := sub< GL(8,k) |
            [ EF_AffineMatrix(rho6(H1.i), [ k | 0 : j in [1..7] ])
              : i in [1..Ngens(H1)] ] cat
            [ EF_AffineMatrix(IdentityMatrix(k,7),
                  EF_Coefficients(b,6)) : b in shifts ] >;
        for v in EF_OrbitRepresentatives(P) do
            if v[8] eq 0 then continue; end if;
            a3 := EF_FormFromVector(R,6,[ v[j] : j in [1..7] ]);
            H := [ <Matrix(g),p> : g in H1, p in [ EF_FormFromVector(R,2,u) : u in V2 ]
                   | EF_WeierstrassChange([ R | a1,0,
                         EF_FormFromVector(R,6,Vector(k,EF_Coefficients(a3,6))*rho6(g)),0,0 ],
                         k!1,p^2+a1*p,p,R!0)[3] eq a3 ];
            Append(~prefixes, <a1,a3,H>);
        end for;
    end for;
    return prefixes;
end function;

function EF2G_Places(R)
    s := R.1;
    t := R.2;
    return [ R | s, t, s+t, s^2+s*t+t^2 ];
end function;

function EF2G_Residues(pi, j, d)
    // Forms of degree d representing every residue class modulo pi^j.
    R := Parent(pi);
    V := VectorSpace(BaseRing(R),d+1);
    W := sub< V | [ V!EF_Coefficients(pi^j*m,d)
                    : m in EF_BinaryBasis(R,d-j*TotalDegree(pi)) ] >;
    return { EF_FormFromVector(R,d,ReduceVector(W,v)) : v in V };
end function;

function EF2G_NonMinimalSubspaces(a1, a3)
    // The (a4,a6) in V8 x V12 for which the model is not minimal at some
    // place pi, as affine subspaces <x,U>: after x -> x+r, y -> y+q*x+v with
    // q, r, v of degrees 2, 4, 6, pi^i divides the new a_i. Only q mod pi,
    // r mod pi^2 and v mod pi^3 matter, and r is fixed modulo pi^2 by the
    // condition on the new a2. For fixed q, r, v the conditions are
    // pi^4 | a4 + c4 and pi^6 | a6 + r*a4 + c6.
    R := Parent(a1);
    V22 := VectorSpace(BaseRing(R),22);
    k := BaseRing(R);
    bad := [];
    for pi in EF2G_Places(R) do
        if not IsDivisibleBy(a1,pi) then continue; end if;
        e := TotalDegree(pi);
        U4 := [ pi^4*m : m in EF_BinaryBasis(R,8-4*e) ];
        U6 := [ pi^6*m : m in EF_BinaryBasis(R,12-6*e) ];
        for q in EF2G_Residues(pi,1,2) do
            r := q*a1+q^2;
            if not IsDivisibleBy(a3+r*a1,pi^3) then continue; end if;
            U := sub< V22 | [ V22!(EF_Coefficients(u,8) cat EF_Coefficients(r*u,12)) : u in U4 ]
                     cat [ V22!(EF_Coefficients(R!0,8) cat EF_Coefficients(u,12)) : u in U6 ] >;
            for v in EF2G_Residues(pi,3,6) do
                c4 := q*a3+(v+r*q)*a1+r^2;
                x := V22!(EF_Coefficients(c4,8) cat EF_Coefficients(r*c4+r^3+v*a3+v^2+r*v*a1,12));
                Append(~bad, <x,U>);
            end for;
        end for;
    end for;
    return bad;
end function;

function FibrationFibreF2G(prefix)
    a1, a3, H := Explode(prefix);
    R := Parent(a1);
    k := BaseRing(R);
    V22 := VectorSpace(k,22);
    G := GL(2,k);
    rho6 := EF_FormAction(R,G,6);
    rho8 := EF_FormAction(R,G,8);
    rho12 := EF_FormAction(R,G,12);
    gens := [];
    for h in H do
        M, p := Explode(h);
        g := G!M;
        q := p^2+a1*p;
        a3M := EF_FormFromVector(R,6,Vector(k,EF_Coefficients(a3,6))*rho6(g));
        B0 := EF_WeierstrassChange([ R | a1,0,a3M,0,0 ],k!1,q,p,R!0);
        rows := [];
        for z in Basis(V22) do
            c := EF_FormFromVector(R,8,Vector(k,[ z[i] : i in [1..9] ])*rho8(g));
            e := EF_FormFromVector(R,12,Vector(k,[ z[i] : i in [10..22] ])*rho12(g));
            B := EF_WeierstrassChange([ R | a1,0,a3M,c,e ],k!1,q,p,R!0);
            Append(~rows, EF_Coefficients(B[4]-B0[4],8) cat EF_Coefficients(B[5]-B0[5],12));
        end for;
        Append(~gens, EF_AffineMatrix(Matrix(k,22,22,&cat rows),
            EF_Coefficients(B0[4],8) cat EF_Coefficients(B0[5],12)));
    end for;
    // The y-shifts y -> y+r form a GF(2)-linear image.
    shifts := [ EF_WeierstrassChange([ R | a1,0,a3,0,0 ],k!1,R!0,R!0,r)
                : r in EF_BinaryBasis(R,6) ];
    S := sub< V22 | [ V22!(EF_Coefficients(B[4],8) cat EF_Coefficients(B[5],12))
                      : B in shifts ] >;
    E, d, toQ, fromQ := EF_AffineQuotient(sub< GL(23,k) | gens >, S);
    W := VectorSpace(k,d);
    // Delta = a3^4 + a1^3*a3^3 + a1^4*a4^2 + a1^5*a3*a4 + a1^6*a6 is affine
    // in (a4,a6) over GF(2).
    linear := Matrix(k,22,25,[ EF_Coefficients(a1^4*c^2+a1^5*a3*c+a1^6*e,24)
        : z in Basis(V22) | true
          where c is EF_FormFromVector(R,8,[ z[i] : i in [1..9] ])
          where e is EF_FormFromVector(R,12,[ z[i] : i in [10..22] ]) ]);
    ok, z0, K := IsConsistent(linear, Vector(k,EF_Coefficients(a3^4+a1^3*a3^3,24)));
    singular := ok select <W!(z0*toQ), sub< W | [ W!(w*toQ) : w in Basis(K) ] >>
                   else false;
    // When a1 = a3 = 0 every model is singular.
    nonminimal := (a1 eq 0 and a3 eq 0) select [] else
        EF_ProjectSubspaces(EF2G_NonMinimalSubspaces(a1,a3), toQ);
    digits := EF_Coefficients(a1,2) cat EF_Coefficients(a3,6);
    return rec< EF_Fibre | Dimension := d, Group := E, Singular := singular,
        NonMinimal := nonminimal,
        Key := func< v | EF_PrimeFieldKey(digits cat Eltseq(Vector(k,Eltseq(v))*fromQ),2) >,
        FromKey := func< key | W!(EF_VectorFromKey(V22, key mod 2^22)*toQ) > >;
end function;
