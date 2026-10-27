// Shared binary-form utilities. Prime-field encodings use most significant
// coefficient first, in the basis s^d,s^(d-1)*t,...,t^d.

function EF_BinaryBasis(R, d)
    assert Rank(R) eq 2 and d ge 0;
    return [ R | R.1^(d-i)*R.2^i : i in [0..d] ];
end function;

function EF_CheckForms(forms, degrees)
    assert #forms eq #degrees and #forms gt 0;
    R := Parent(forms[1]);
    assert Type(R) eq RngMPol and Rank(R) eq 2;
    assert Type(BaseRing(R)) eq FldFin;
    for i in [1..#forms] do
        assert Parent(forms[i]) eq R;
        assert &and[ TotalDegree(m) eq degrees[i] : m in Monomials(forms[i]) ];
    end for;
    return R;
end function;

function EF_Coefficients(f, d)
    return [ MonomialCoefficient(f,m) : m in EF_BinaryBasis(Parent(f),d) ];
end function;

function EF_FormFromVector(R, d, v)
    assert #Eltseq(v) eq d+1;
    basis := EF_BinaryBasis(R,d);
    return &+[ R | v[i]*basis[i] : i in [1..d+1] ];
end function;

function EF_BaseChange(f, M)
    R := Parent(f);
    return Evaluate(f,[M[1,1]*R.1+M[1,2]*R.2,
                       M[2,1]*R.1+M[2,2]*R.2]);
end function;

function EF_WeierstrassChange(A, u, r, s, t)
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

function EF_PrimeFieldKey(v, p)
    assert IsPrime(p);
    n := 0;
    for a in Eltseq(v) do
        digit := Integers()!a;
        assert digit ge 0 and digit lt p;
        n := p*n+digit;
    end for;
    return n;
end function;

function EF_PrimeFieldKeys(M, p)
    weights := Matrix(Integers(),Ncols(M),1,
                      [ p^(Ncols(M)-i) : i in [1..Ncols(M)] ]);
    return Eltseq(ChangeRing(M,Integers())*weights);
end function;

function EF_VectorFromKey(V, n)
    p := #BaseRing(V);
    assert IsPrime(p);
    d := Dimension(V);
    assert n ge 0 and n lt p^d;
    return V![ (n div p^(d-i)) mod p : i in [1..d] ];
end function;

function EF_FormFromKey(R, d, n)
    return EF_FormFromVector(R,d,
        EF_VectorFromKey(VectorSpace(BaseRing(R),d+1),n));
end function;

function EF_ReducedComplement(V, S)
    // Image of the canonical coset reduction. No invariant complement needed.
    Q := sub< V | [ ReduceVector(S,V.i) : i in [1..Dimension(V)] ] >;
    assert Dimension(Q)+Dimension(S) eq Dimension(V);
    assert Dimension(Q meet S) eq 0;
    return Q;
end function;

function EF_AffineFixedCount(L, b)
    // Row-vector action v -> v*L+b. Fixed points solve v*(I-L)=b.
    k := BaseRing(L);
    d := Nrows(L);
    assert Ncols(L) eq d and #Eltseq(b) eq d;
    D := IdentityMatrix(k,d)-L;
    r := Rank(D);
    if Rank(VerticalJoin(D,Matrix(k,1,d,Eltseq(b)))) ne r then
        return 0;
    end if;
    return (#k)^(d-r);
end function;
