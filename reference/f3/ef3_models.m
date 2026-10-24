// Completed-square models are triples [C,A,B], not long Weierstrass lists.
// Load ef_common.m first.

function EF3_CheckModel(model)
    R := EF_CheckForms(model,[4,8,12]);
    assert #BaseRing(R) eq 3;
    return R;
end function;

function EF3_Key(v)
    return EF_PrimeFieldKey(v,3);
end function;

function EF3_VectorFromKey(V, n)
    assert #BaseRing(V) eq 3;
    return EF_VectorFromKey(V,n);
end function;

function EF3_FormFromKey(R, d, n)
    assert #BaseRing(R) eq 3;
    return EF_FormFromKey(R,d,n);
end function;

function EF3_ModelKey(model)
    R := EF3_CheckModel(model);
    return EF3_Key(EF_Coefficients(model[1],4) cat
                   EF_Coefficients(model[2],8) cat
                   EF_Coefficients(model[3],12));
end function;

function FibrationFromKeyF3(R, key)
    assert #BaseRing(R) eq 3;
    assert key ge 0 and key lt 3^27;
    C := EF3_FormFromKey(R,4,key div 3^22);
    A := EF3_FormFromKey(R,8,(key div 3^13) mod 3^9);
    B := EF3_FormFromKey(R,12,key mod 3^13);
    return [ R | C,A,B ];
end function;

function EF3_AsLongModel(model)
    R := EF3_CheckModel(model);
    return [ R | 0,model[1],0,model[2],model[3] ];
end function;

function EF3_NormalizeLongModel(model)
    R := EF_CheckForms(model,[2,4,6,8,12]);
    assert #BaseRing(R) eq 3;
    a1,a2,a3,a4,a6 := Explode(model);
    return [ R | a2+a1^2/4, a4+a1*a3/2, a6+a3^2/4 ];
end function;

function FibrationDiscriminantF3(model)
    R := EF3_CheckModel(model);
    C,A,B := Explode(model);
    return C^2*A^2-A^3-C^3*B;
end function;

function EF3_Transform(model, M, q : u := 1)
    // Base substitution, then scaling, then x -> x+q.
    R := EF3_CheckModel(model);
    k := BaseRing(R);
    assert BaseRing(M) eq k and Nrows(M) eq 2 and Ncols(M) eq 2;
    assert Determinant(M) ne 0;
    u := k!u;
    assert u ne 0;
    q := R!q;
    dummy := EF_CheckForms([q],[4]);
    C := u^2*EF_BaseChange(model[1],M);
    A := u^4*EF_BaseChange(model[2],M);
    B := u^6*EF_BaseChange(model[3],M);
    return [ R | C,A-C*q,B+A*q+C*q^2+q^3 ];
end function;

function EF3_CanonicalBaseMatrix(M)
    for a in Eltseq(M) do
        if a ne 0 then return a^-1*M; end if;
    end for;
    error "The zero matrix has no projective representative";
end function;

function EF3_BaseMatrices(k)
    assert #k eq 3;
    V := VectorSpace(k,4);
    matrices := [ MatrixAlgebra(k,2) | ];
    for n in [0..80] do
        M := Matrix(k,2,2,Eltseq(EF3_VectorFromKey(V,n)));
        if Determinant(M) ne 0 and EF3_CanonicalBaseMatrix(M) eq M then
            Append(~matrices,M);
        end if;
    end for;
    assert #matrices eq 24;
    return matrices;
end function;

function EF3_Compose(g, h)
    // Effective elements are pairs <M,q>. Apply g, then h.
    // Pullbacks satisfy (f(M*z))(N*z)=f(M*N*z).
    M,q := Explode(g);
    N,r := Explode(h);
    return <EF3_CanonicalBaseMatrix(M*N),EF_BaseChange(q,N)+r>;
end function;

function EF3_Inverse(g)
    M,q := Explode(g);
    N := M^-1;
    return <EF3_CanonicalBaseMatrix(N),-EF_BaseChange(q,N)>;
end function;
