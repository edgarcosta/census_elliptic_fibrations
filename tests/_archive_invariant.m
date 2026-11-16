function form_vector(f, d)
    R := Parent(f);
    return [BaseRing(R) | MonomialCoefficient(f, R.1^(d-j)*R.2^j) : j in [0..d]];
end function;

function vector_form(R, d, v)
    return &+[R | v[j+1]*R.1^(d-j)*R.2^j : j in [0..d]];
end function;

function number(v, q)
    n := 0;
    for c in Eltseq(v) do n := q*n+Integers()!c; end for;
    return n;
end function;

function prefix_vector(W)
    return Characteristic(Universe(W)) eq 2
        select form_vector(W[1], 2) cat form_vector(W[3], 6)
        else form_vector(W[2], 4) cat form_vector(W[4], 8);
end function;

function suffix_vector(W)
    return Characteristic(Universe(W)) eq 2
        select form_vector(W[4], 8) cat form_vector(W[5], 12)
        else form_vector(W[5], 12);
end function;

function with_suffix(W, v)
    R := Universe(W);
    if Characteristic(R) eq 2 then
        return [R | W[1], 0, W[3], vector_form(R, 8, v[1..9]),
                    vector_form(R, 12, v[10..22])];
    end if;
    return [R | 0, W[2], 0, W[4], vector_form(R, 12, v)];
end function;

function change(W, M, h, z)
    R := Universe(W);
    a1, a2, a3, a4, a6 := Explode([R | Evaluate(f,
        [M[1, 1]*R.1+M[1, 2]*R.2, M[2, 1]*R.1+M[2, 2]*R.2]) : f in W]);
    if Characteristic(R) eq 2 then
        r := h^2+a1*h;
        return [R | a1, 0, a3+r*a1,
            a4-h*a3-(z+r*h)*a1+3*r^2-2*h*z,
            a6+r*a4+r^3-z*a3-r*z*a1-z^2];
    end if;
    return [R | 0, a2, 0, a4-a2*h, a6+a4*h+a2*h^2+h^3];
end function;

CONTEXT := recformat<Prefix, Width, Matrices, Offsets>;
function invariant_context(W)
    R := Universe(W);
    q := Characteristic(R);
    d := q eq 2 select 22 else 13;
    V := VectorSpace(GF(q), d);
    W := with_suffix(W, Eltseq(V!0));
    transports := [* *];
    best := q^(q eq 2 select 10 else 14);
    for g in GL(2, GF(q)) do
        M := Matrix(g);
        if q eq 3 then
            first := [c : c in Eltseq(M) | c ne 0][1];
            if first ne 1 then continue; end if;
        end if;
        if q eq 2 then
            hs := [vector_form(R, 2, Eltseq(v)) : v in VectorSpace(GF(q), 3)];
        else
            T := change(W, M, R!0, R!0);
            h := R!0;
            if T[2] ne 0 then
                U := sub<VectorSpace(GF(q), 9) |
                    [form_vector(T[2]*R.1^(4-j)*R.2^j, 8) : j in [0..4]]>;
                A := vector_form(R, 8, Eltseq(ReduceVector(U,
                    Vector(GF(q), form_vector(T[4], 8)))));
                ok, h := IsDivisibleBy(T[4]-A, T[2]);
                assert ok;
            end if;
            hs := [h];
        end if;
        for h in hs do
            Z := change(W, M, h, R!0);
            n := number(prefix_vector(Z), q);
            if n lt best then
                best := n;
                target := Z;
                transports := [* *];
            end if;
            if n eq best then Append(~transports, <M, h>); end if;
        end for;
    end for;
    if q eq 2 then
        U := sub<V | [suffix_vector(change(with_suffix(target, Eltseq(V!0)),
            IdentityMatrix(GF(q), 2), R!0, R.1^(6-j)*R.2^j)) : j in [0..6]]>;
    elif target[2] eq 0 then
        U := sub<V | [suffix_vector(change(with_suffix(target, Eltseq(V!0)),
            IdentityMatrix(GF(q), 2), R.1^(4-j)*R.2^j, R!0)) : j in [0..4]]>;
    else
        U := sub<V |>;
    end if;
    matrices := [* *];
    offsets := [* *];
    for transport in transports do
        M, h := Explode(transport);
        offset := ReduceVector(U, V!suffix_vector(change(W, M, h, R!0)));
        rows := [ReduceVector(U, V!suffix_vector(change(with_suffix(W, Eltseq(v)),
            M, h, R!0)))-offset : v in Basis(V)];
        L := Matrix(GF(q), d, d, &cat[Eltseq(v) : v in rows]);
        if &or[L eq matrices[j] and offset eq offsets[j] : j in [1..#matrices]] then
            continue;
        end if;
        Append(~matrices, L);
        Append(~offsets, offset);
    end for;
    return rec<CONTEXT | Prefix := best, Width := d, Matrices := matrices,
        Offsets := offsets>;
end function;

function orbit_invariant(W, C)
    q := Characteristic(Universe(W));
    v := Vector(GF(q), suffix_vector(W));
    tail := Min([number(v*C`Matrices[j]+C`Offsets[j], q) : j in [1..#C`Matrices]]);
    return C`Prefix*q^C`Width+tail;
end function;
