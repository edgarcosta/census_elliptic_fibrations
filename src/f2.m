import "forms.m": binary_basis, coefficients, form_from_vector, form_action;
import "models.m": weierstrass_change, prefix_key_map;
import "engine.m": orbit_representatives, affine_matrix, orbit_fibre;
import "minimality.m": nonminimal_locus;
import "classification.m": generic_loci;

function normalized_shift(a1, q)
    return q^2+a1*q;
end function;

function f2_prefixes(R)
    k := BaseRing(R);
    assert #k eq 2 and Rank(R) eq 2;
    G := GL(2,k);
    rho2 := form_action(R,G,2);
    rho6 := form_action(R,G,6);
    G2 := Image(rho2);
    B2 := binary_basis(R,2);
    V2 := VectorSpace(k,3);
    prefixes := [* *];
    for c in orbit_representatives(G2) do
        a1 := form_from_vector(R,2,c);
        H1 := Stabiliser(G2,c) @@ rho2;
            shifts := [ weierstrass_change([ R | a1,0,0,0,0 ],
                        k!1,normalized_shift(a1,m),m,R!0)[3] : m in B2 ];
        P := sub< GL(8,k) |
            [ affine_matrix(rho6(H1.i), [ k | 0 : j in [1..7] ])
              : i in [1..Ngens(H1)] ] cat
            [ affine_matrix(IdentityMatrix(k,7),
                  coefficients(b,6)) : b in shifts ] >;
        for v in orbit_representatives(P) do
            if v[8] eq 0 then continue; end if;
            a3 := form_from_vector(R,6,[ v[j] : j in [1..7] ]);
            H := [ <Matrix(g),p> : g in H1, p in [ form_from_vector(R,2,u) : u in V2 ]
                   | weierstrass_change([ R | a1,0,
                         form_from_vector(R,6,Vector(k,coefficients(a3,6))*rho6(g)),0,0 ],
                         k!1,normalized_shift(a1,p),p,R!0)[3] eq a3 ];
            Append(~prefixes, <a1,a3,H>);
        end for;
    end for;
    return prefixes;
end function;

function f2_fibre(prefix)
    a1, a3, H := Explode(prefix);
    R := Parent(a1);
    k := BaseRing(R);
    V22 := VectorSpace(k,22);
    G := GL(2,k);
    rho6 := form_action(R,G,6);
    rho8 := form_action(R,G,8);
    rho12 := form_action(R,G,12);
    gens := [];
    for h in H do
        M, p := Explode(h);
        g := G!M;
        q := normalized_shift(a1,p);
        a3M := form_from_vector(R,6,Vector(k,coefficients(a3,6))*rho6(g));
        B0 := weierstrass_change([ R | a1,0,a3M,0,0 ],k!1,q,p,R!0);
        rows := [];
        for z in Basis(V22) do
            c := form_from_vector(R,8,Vector(k,[ z[i] : i in [1..9] ])*rho8(g));
            e := form_from_vector(R,12,Vector(k,[ z[i] : i in [10..22] ])*rho12(g));
            B := weierstrass_change([ R | a1,0,a3M,c,e ],k!1,q,p,R!0);
            Append(~rows, coefficients(B[4]-B0[4],8) cat coefficients(B[5]-B0[5],12));
        end for;
        Append(~gens, affine_matrix(Matrix(k,22,22,&cat rows),
            coefficients(B0[4],8) cat coefficients(B0[5],12)));
    end for;
    shifts := [ weierstrass_change([ R | a1,0,a3,0,0 ],k!1,R!0,R!0,r)
                : r in binary_basis(R,6) ];
    S := sub< V22 | [ V22!(coefficients(B[4],8) cat coefficients(B[5],12))
                      : B in shifts ] >;
    S0,R0 := generic_loci(R,[a1,a3]);
    N := nonminimal_locus(R,[a1,a3]);
    return orbit_fibre(V22,sub< GL(23,k) | gens >,S,[* S0,R0,N *],
        prefix_key_map(R,[a1,a3]));
end function;
