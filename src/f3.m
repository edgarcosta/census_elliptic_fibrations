import "forms.m": binary_basis, coefficients, form_from_vector, form_action;
import "models.m": weierstrass_change, prefix_key_map;
import "engine.m": orbit_representatives, affine_matrix, orbit_fibre;
import "minimality.m": nonminimal_locus;
import "classification.m": generic_loci;

function f3_prefixes(R)
    k := BaseRing(R);
    assert #k eq 3 and Rank(R) eq 2;
    G := GL(2,k);
    rho4 := form_action(R,G,4);
    rho8 := form_action(R,G,8);
    G4 := Image(rho4);
    G8 := Image(rho8);
    B4 := binary_basis(R,4);
    prefixes := [* *];
    for c in orbit_representatives(G4) do
        C := form_from_vector(R,4,c);
        if C eq 0 then
            for a in orbit_representatives(G8) do
                H := Stabiliser(G8,a) @@ rho8;
                Append(~prefixes, <C,form_from_vector(R,8,a),H>);
            end for;
            continue;
        end if;
        HC := Stabiliser(G4,c) @@ rho4;
        M8 := GModule(HC, [ rho8(HC.i) : i in [1..Ngens(HC)] ]);
        shifts := [ weierstrass_change([ R | 0,C,0,0,0 ],k!1,m,R!0,R!0)[4]
                    : m in B4 ];
        S := sub< M8 | [ M8!coefficients(b,8) : b in shifts ] >;
        assert Dimension(S) eq 5;
        Q, proj := quo< M8 | S >;
        phi := hom< HC -> GL(4,k) |
            [ GL(4,k)!X : X in ActionGenerators(Q) ] >;
        GQ := Image(phi);
        for abar in orbit_representatives(GQ) do
            A := form_from_vector(R,8,Eltseq((Q!abar) @@ proj));
            Append(~prefixes, <C,A,Stabiliser(GQ,abar) @@ phi>);
        end for;
    end for;
    return prefixes;
end function;

function f3_fibre(prefix)
    C, A, H := Explode(prefix);
    R := Parent(C);
    k := BaseRing(R);
    V12 := VectorSpace(k,13);
    rho8 := form_action(R,H,8);
    rho12 := form_action(R,H,12);
    if C ne 0 then
        gens := [];
        for i in [1..Ngens(H)] do
            AM := form_from_vector(R,8,Vector(k,coefficients(A,8))*rho8(H.i));
            ok, q := IsDivisibleBy(AM-A, C);
            assert ok;
            B := weierstrass_change([ R | 0,C,0,AM,0 ],k!1,q,R!0,R!0);
            Append(~gens, affine_matrix(rho12(H.i),
                coefficients(B[5],12)));
        end for;
        E := sub< GL(14,k) | gens >;
        assert 2*#E eq #H;
        S := sub< V12 | >;
    else
        E12 := sub< GL(14,k) | [ affine_matrix(rho12(H.i), [ k | 0 : j in [1..13] ])
                                 : i in [1..Ngens(H)] ] >;
        images := [ coefficients(weierstrass_change(
                        [ R | 0,0,0,A,0 ],k!1,m,R!0,R!0)[5],12)
                    : m in binary_basis(R,4) ];
        S := sub< V12 | [ V12!w : w in images ] >;
        E := E12;
    end if;
    S0,R0 := generic_loci(R,[C,A]);
    N := nonminimal_locus(R,[C,A]);
    return orbit_fibre(V12,E,S,[* S0,R0,N *],prefix_key_map(R,[C,A]));
end function;
