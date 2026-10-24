function binary_basis(R, d)
    assert Rank(R) eq 2;
    return [ R | R.1^(d-i)*R.2^i : i in [0..d] ];
end function;

function coefficients(f, d)
    return [ MonomialCoefficient(f,m) : m in binary_basis(Parent(f),d) ];
end function;

function form_from_vector(R, d, v)
    assert #Eltseq(v) eq d+1;
    basis := binary_basis(R,d);
    return &+[ R | v[i]*basis[i] : i in [1..d+1] ];
end function;

function base_change(f, M)
    R := Parent(f);
    return Evaluate(f,[M[1,1]*R.1+M[1,2]*R.2,
                       M[2,1]*R.1+M[2,2]*R.2]);
end function;

function form_action(R, G, d)
    assert BaseRing(G) eq BaseRing(R);
    V := SymmetricPower(GModule(G),d);
    return hom< G -> GL(d+1,BaseRing(R)) |
        ActionGenerators(V) >;
end function;

function is_form_power(f, n)
    return &and[ &and[ e mod n eq 0 : e in Exponents(m) ] : m in Monomials(f) ];
end function;
