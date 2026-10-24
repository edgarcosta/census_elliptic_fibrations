import "forms.m": binary_basis, coefficients, form_from_vector;

function weierstrass_change(A, u, r, s, t)
    // Substitute x_old = u^2*x+r, y_old = u^3*y+s*u^2*x+t.
    // These formulas also apply when the discriminant is zero.
    a1, a2, a3, a4, a6 := Explode(A);
    v := u^-1;
    return [ Universe(A) |
        v*(a1+2*s),
        v^2*(a2-s*a1+3*r-s^2),
        v^3*(a3+r*a1+2*t),
        v^4*(a4-s*a3+2*r*a2-(t+r*s)*a1+3*r^2-2*s*t),
        v^6*(a6+r*a4+r^2*a2+r^3-t*a3-r*t*a1-t^2) ];
end function;

function prime_field_key(v, p)
    assert IsPrime(p);
    n := 0;
    for a in Eltseq(v) do
        digit := Integers()!a;
        assert digit ge 0 and digit lt p;
        n := p*n+digit;
    end for;
    return n;
end function;

function prime_field_keys(M, p)
    weights := Matrix(Integers(),Ncols(M),1,
                      [ p^(Ncols(M)-i) : i in [1..Ncols(M)] ]);
    return Eltseq(ChangeRing(M,Integers())*weights);
end function;

function vector_from_key(V, n)
    p := #BaseRing(V);
    assert IsPrime(p);
    d := Dimension(V);
    assert n ge 0 and n lt p^d;
    return V![ (n div p^(d-i)) mod p : i in [1..d] ];
end function;

function form_from_key(R, d, n)
    return form_from_vector(R,d,
        vector_from_key(VectorSpace(BaseRing(R),d+1),n));
end function;

function valid_ring(R)
    return Type(R) eq RngMPol and Rank(R) eq 2
        and Type(BaseRing(R)) eq FldFin and #BaseRing(R) in {2,3};
end function;

function valid_model(W)
    if #W ne 5 or not valid_ring(Universe(W)) then return false; end if;
    R := Universe(W);
    degrees := [2,4,6,8,12];
    return &and[ Parent(W[i]) eq R and
        &and[ TotalDegree(m) eq degrees[i] : m in Monomials(W[i]) ]
        : i in [1..5] ];
end function;

function checked_model(W)
    assert valid_model(W);
    return Universe(W);
end function;

function normalize_model(W)
    R := Universe(W);
    a1,a2,a3,a4,a6 := Explode(W);
    if Characteristic(R) eq 2 then
        return [ R | a1,0,a3+a1*a2,a4+a2^2,a6+a2*a4 ];
    end if;
    return [ R | 0,a2+a1^2/4,0,a4+a1*a3/2,a6+a3^2/4 ];
end function;

function model_discriminant(W)
    a1,a2,a3,a4,a6 := Explode(normalize_model(W));
    if Characteristic(Universe(W)) eq 2 then
        return a3^4+a1^3*a3^3+a1^4*a4^2+a1^5*a3*a4+a1^6*a6;
    end if;
    return a2^2*a4^2-a4^3-a2^3*a6;
end function;

function model_key(W)
    W := normalize_model(W);
    indices := Characteristic(Universe(W)) eq 2 select [1,3,4,5] else [2,4,5];
    degrees := [2,4,6,8,12];
    return prime_field_key(&cat[ coefficients(W[i],degrees[i]) : i in indices ],
        Characteristic(Universe(W)));
end function;

function model_from_key(R, key)
    p := Characteristic(R);
    indices := p eq 2 select [1,3,4,5] else [2,4,5];
    degrees := [2,4,6,8,12];
    W := [ R | 0,0,0,0,0 ];
    for i in Reverse(indices) do
        W[i] := form_from_key(R,degrees[i],key mod p^(degrees[i]+1));
        key div:= p^(degrees[i]+1);
    end for;
    return W;
end function;

function prefix_key_map(R, forms)
    p := Characteristic(R);
    degrees := p eq 2 select [2,6] else [4,8];
    width := p eq 2 select 22 else 13;
    digits := &cat[ coefficients(forms[i],degrees[i]) : i in [1,2] ];
    prefix := prime_field_key(digits,p)*p^width;
    return func< P | [ Integers() | prefix+b : b in prime_field_keys(P,p) ] >;
end function;

intrinsic aInvariantsFromKey(R::RngMPol, key::RngIntElt) -> SeqEnum
{The normalized model [a1,a2,a3,a4,a6] over R encoded by key, with
0 <= key < 2^32 over GF(2) or 0 <= key < 3^27 over GF(3).
See also FibrationKey, FibrationProperties}
    require valid_ring(R): "R must be a binary polynomial ring over GF(2) or GF(3)";
    bound := Characteristic(R) eq 2 select 2^32 else 3^27;
    require 0 le key and key lt bound: "key is outside the encoding range";
    return model_from_key(R,key);
end intrinsic;

intrinsic FibrationKey(W::SeqEnum) -> RngIntElt
{The coefficient key of the normalized five-coefficient Weierstrass model W.
Normalization removes a2 over GF(2) and completes the square over GF(3).
See also aInvariantsFromKey}
    require valid_model(W): "W must have five homogeneous coefficients of degrees 2,4,6,8,12 over one supported ring";
    return model_key(W);
end intrinsic;
