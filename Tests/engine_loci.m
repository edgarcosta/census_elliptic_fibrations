import "src/engine.m": orbit_fibre, affine_matrix, fibre_counts, fibre_representatives;

// Arbitrary affine groups exceed the public census domain; brute force checks
// the engine's invariant unions and shift quotient without using its counters.
procedure check_engine(V, E, T, loci)
    p := #BaseRing(V);
    d := Dimension(V);
    keys := func< P | [ Integers() | &+[ Integers()!P[i,j]*p^(d-j) : j in [1..d] ] : i in [1..Nrows(P)] ] >;
    F := orbit_fibre(V,E,T,loci,keys);
    expected := AssociativeArray();
    for flags in CartesianPower([false,true],#loci) do expected[<x : x in flags>] := 0; end for;
    seen := {};
    orbit_of := AssociativeArray();
    next := 0;
    for v in V do
        if v in seen then continue; end if;
        next +:= 1;
        orbit := { V![ (Vector(BaseRing(V),Eltseq(v) cat [1])*Matrix(g))[j] : j in [1..d] ]+t
            : g in E, t in T };
        pattern := < &or[ (v-A[1]) in A[2] : A in union ] : union in loci >;
        expected[pattern] +:= 1;
        for w in orbit do
            number := &+[ Integers()!w[j]*p^(d-j) : j in [1..d] ];
            orbit_of[number] := next;
        end for;
        seen join:= orbit;
    end for;
    actual := fibre_counts(F);
    assert Keys(actual) eq Keys(expected);
    for pattern in Keys(expected) do assert actual[pattern] eq expected[pattern]; end for;
    representatives := fibre_representatives(F,Setseq(Keys(expected)));
    assert #representatives eq next and #representatives gt 0;
    assert #{ orbit_of[pair[1]] : pair in representatives } eq next;
    for pair in representatives do
        v := V![ (pair[1] div p^(d-j)) mod p : j in [1..d] ];
        assert pair[2] eq < &or[ (v-A[1]) in A[2] : A in union ] : union in loci >;
    end for;
end procedure;

// Pins empty/full loci, overlapping components and invariant unions whose
// components are exchanged by translations; quotient lifts retain orbit identity.
procedure test_affine_loci()
    V := VectorSpace(GF(2),3);
    E := sub< GL(4,GF(2)) | affine_matrix(IdentityMatrix(GF(2),3),V.3) >;
    T := sub< V | V.2 >;
    check_engine(V,E,T,[* *]);
    check_engine(V,E,V,[* [* <V!0,V> *] *]);
    check_engine(V,E,T,[* [* *],[* <V!0,V> *] *]);
    U := sub< V | V.2,V.3 >;
    A := sub< V | V.2 >;
    union := [* <V!0,A>,<V.3,A>,<V!0,U> *];
    check_engine(V,E,T,[* union,[* <V!0,V> *],union *]);
end procedure;

// Pins batch enumeration for an affine action without a fixed point.
procedure test_large_affine_listing()
    k := GF(3);
    d := 11;
    V := VectorSpace(k,d);
    E := sub< GL(d+1,k) | affine_matrix(IdentityMatrix(k,d),V.d) >;
    loci := [* [* <V!0,sub< V | [ V.i : i in [2..d] ] >> *],
        [* <V.2,sub< V | [ V.i : i in [1..d] | i ne 2 ] >> *] *];
    keys := func< P | [ Integers() | &+[ Integers()!P[i,j]*3^(d-j) : j in [1..d] ] : i in [1..Nrows(P)] ] >;
    F := orbit_fibre(V,E,sub< V | >,loci,keys);
    reps := fibre_representatives(F,[<false,false>]);
    assert #reps eq 4*3^8;
    assert { pair[1] : pair in reps } eq { n : n in [0..3^d-1] |
        n div 3^(d-1) ne 0 and (n div 3^(d-2)) mod 3 ne 1 and n mod 3 eq 0 };
end procedure;

test_affine_loci();
test_large_affine_listing();

// Pins overlap correction for large affine spaces without listing their points.
procedure test_large_union_count()
    k := GF(2);
    V := VectorSpace(k,20);
    E := sub< GL(21,k) | >;
    U := sub< V | [ V.i : i in [2..20] ] >;
    W := sub< V | [ V.i : i in [1..20] | i ne 2 ] >;
    F := orbit_fibre(V,E,sub< V | >,[* [* <V!0,U>,<V!0,W> *] *],
        func< P | [ Integers() | 0 : i in [1..Nrows(P)] ] >);
    counts := fibre_counts(F);
    assert counts[<true>] eq 3*2^18 and counts[<false>] eq 2^18;
end procedure;

test_large_union_count();
