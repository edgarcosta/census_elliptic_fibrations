// The census over GF(3) for ef_orbits.m: models y^2 = x^3 + C*x^2 + A*x + B
// with [C,A,B] of degrees 4, 8, 12. Load ef_common.m, ef3_models.m and
// ef_orbits.m first. A prefix is <C,A,H> with H <= GL(2,3) the stabilizer
// of C and of the class of A modulo C*V4 (C != 0), or of A (C = 0).

function FibrationOrbitPrefixesF3G(R)
    k := BaseRing(R);
    assert #k eq 3 and Rank(R) eq 2;
    G := GL(2,k);
    rho4 := EF_FormAction(R,G,4);
    rho8 := EF_FormAction(R,G,8);
    G4 := Image(rho4);
    G8 := Image(rho8);
    B4 := EF_BinaryBasis(R,4);
    prefixes := [* *];
    for c in EF_OrbitRepresentatives(G4) do
        C := EF_FormFromVector(R,4,c);
        if C eq 0 then
            for a in EF_OrbitRepresentatives(G8) do
                H := Stabiliser(G8,a) @@ rho8;
                Append(~prefixes, <C,EF_FormFromVector(R,8,a),H>);
            end for;
            continue;
        end if;
        HC := Stabiliser(G4,c) @@ rho4;
        M8 := GModule(HC, [ rho8(HC.i) : i in [1..Ngens(HC)] ]);
        shifts := [ EF_WeierstrassChange([ R | 0,C,0,0,0 ],k!1,m,R!0,R!0)[4]
                    : m in B4 ];
        S := sub< M8 | [ M8!EF_Coefficients(b,8) : b in shifts ] >;
        assert Dimension(S) eq 5;
        Q, proj := quo< M8 | S >;
        phi := hom< HC -> GL(4,k) |
            [ GL(4,k)!X : X in ActionGenerators(Q) ] >;
        GQ := Image(phi);
        for abar in EF_OrbitRepresentatives(GQ) do
            A := EF_FormFromVector(R,8,Eltseq((Q!abar) @@ proj));
            Append(~prefixes, <C,A,Stabiliser(GQ,abar) @@ phi>);
        end for;
    end for;
    return prefixes;
end function;

function EF3G_Places(R)
    s := R.1;
    t := R.2;
    return [ R | s, t, s+t, s-t, s^2+t^2, s^2+s*t-t^2, s^2-s*t-t^2 ];
end function;

function EF3G_NonMinimalSubspaces(C, A)
    // The B in V12 for which [C,A,B] is not minimal at some place pi, as
    // cosets <x,U>: pi^2 | C, and some r in V4 gives pi^4 | A-C*r and
    // pi^6 | B+A*r+C*r^2+r^3, a coset of U = pi^6*V_(12-6*deg pi).
    R := Parent(C);
    k := BaseRing(R);
    V4 := VectorSpace(k,5);
    V12 := VectorSpace(k,13);
    bad := [];
    for pi in EF3G_Places(R) do
        e := TotalDegree(pi);
        if not IsDivisibleBy(C,pi^2) then continue; end if;
        if C eq 0 and not IsDivisibleBy(A,pi^4) then continue; end if;
        P6 := sub< V12 | [ V12!EF_Coefficients(pi^6*m,12)
                           : m in EF_BinaryBasis(R,12-6*e) ] >;
        cosets := {@ V12 | @};
        for v in V4 do
            r := EF_FormFromVector(R,4,v);
            if IsDivisibleBy(A-C*r,pi^4) then
                Include(~cosets, ReduceVector(P6,
                    V12!EF_Coefficients(-(A*r+C*r^2+r^3),12)));
            end if;
        end for;
        bad cat:= [ <c,P6> : c in cosets ];
    end for;
    return bad;
end function;

function FibrationFibreF3G(prefix)
    C, A, H := Explode(prefix);
    R := Parent(C);
    k := BaseRing(R);
    V12 := VectorSpace(k,13);
    rho8 := EF_FormAction(R,H,8);
    rho12 := EF_FormAction(R,H,12);
    if C ne 0 then
        gens := [];
        for i in [1..Ngens(H)] do
            AM := EF_FormFromVector(R,8,Vector(k,EF_Coefficients(A,8))*rho8(H.i));
            // The x-translation that restores A after M.
            ok, q := IsDivisibleBy(AM-A, C);
            assert ok;
            B := EF_WeierstrassChange([ R | 0,C,0,AM,0 ],k!1,q,R!0,R!0);
            Append(~gens, EF_AffineMatrix(rho12(H.i),
                EF_Coefficients(B[5],12)));
        end for;
        E := sub< GL(14,k) | gens >;
        // Only -I acts trivially, since q is determined by M.
        assert 2*#E eq #H;
        d := 13;
        toQ := IdentityMatrix(k,13);
        fromQ := IdentityMatrix(k,13);
        // Delta = C^2*A^2 - A^3 - C^3*B is affine in B with a unique root.
        ok, B0 := IsDivisibleBy(A^2*(C^2-A), C^3);
        singular := ok select <V12!EF_Coefficients(B0,12), sub< V12 | >>
                       else false;
    else
        E12 := sub< GL(14,k) | [ EF_AffineMatrix(rho12(H.i), [ k | 0 : j in [1..13] ])
                                 : i in [1..Ngens(H)] ] >;
        images := [ EF_Coefficients(EF_WeierstrassChange(
                        [ R | 0,0,0,A,0 ],k!1,m,R!0,R!0)[5],12)
                    : m in EF_BinaryBasis(R,4) ];
        S := sub< V12 | [ V12!w : w in images ] >;
        E, d, toQ, fromQ := EF_AffineQuotient(E12, S);
        W := VectorSpace(k,d);
        // Delta = -A^3.
        singular := A eq 0 select <W!0, W> else false;
    end if;
    W := VectorSpace(k,d);
    nonminimal := (C eq 0 and A eq 0) select [] else
        EF_ProjectSubspaces(EF3G_NonMinimalSubspaces(C,A), toQ);
    digits := EF_Coefficients(C,4) cat EF_Coefficients(A,8);
    prefixKey := EF3_Key(digits)*3^13;
    return rec< EF_Fibre | Dimension := d, Group := E, Singular := singular,
        NonMinimal := nonminimal,
        Key := func< v | EF3_Key(digits cat Eltseq(Vector(k,Eltseq(v))*fromQ)) >,
        Keys := func< P | [ prefixKey+b : b in EF_PrimeFieldKeys(P*fromQ,3) ] >,
        FromKey := func< key | W!(Vector(k,EF_Coefficients(
            FibrationFromKeyF3(R,key)[3],12))*toQ) > >;
end function;
