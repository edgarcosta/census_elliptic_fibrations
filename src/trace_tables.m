TRACE_TABLE_FORMAT := recformat<Polynomial, Field, Generator, Ordinary, NodalTrace>;

function trace_polynomial(p, d)
    polys2 := [
        [1, 1], [1, 1, 1], [1, 1, 0, 1], [1, 1, 0, 0, 1], [1, 0, 1, 0, 0, 1],
        [1, 1, 0, 1, 1, 0, 1], [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 1, 1, 1, 0, 0, 0, 1],
        [1, 0, 0, 0, 1, 0, 0, 0, 0, 1], [1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1]];
    polys3 := [
        [1, 1], [2, 2, 1], [1, 2, 0, 1], [2, 0, 0, 2, 1], [1, 2, 0, 0, 0, 1],
        [2, 2, 1, 0, 2, 0, 1], [1, 0, 2, 0, 0, 0, 0, 1], [2, 2, 2, 0, 1, 2, 0, 0, 1],
        [1, 1, 2, 2, 0, 0, 0, 0, 0, 1], [2, 1, 0, 0, 2, 2, 2, 0, 0, 0, 1],
        [1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1], [2, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1]];
    P := PolynomialRing(GF(p));
    return P!(p eq 2 select polys2[d] else polys3[d]);
end function;

function coordinate_index(a, weights)
    digits := Eltseq(a);
    return 1+&+[Integers()!digits[i]*weights[i] : i in [1..#digits]];
end function;

// Coordinates form (Z/3Z)^d; no twiddle factors or digit reversal are needed.
procedure radix3(~v, w)
    w2 := w^2;
    stride := 1;
    while stride lt #v do
        for block in [0..(#v div (3*stride))-1], j in [0..stride-1] do
            pos := block*3*stride+j+1;
            x := v[pos]; y := v[pos+stride]; z := v[pos+2*stride];
            v[pos] := x+y+z;
            v[pos+stride] := x+w*y+w2*z;
            v[pos+2*stride] := x+w2*y+w*z;
        end for;
        stride *:= 3;
    end while;
end procedure;

function ordinary3(K, alpha, d)
    Q := #K;
    M := 1062907;
    R := Integers(M);
    w := R!662797;
    assert IsPrime(M) and M mod 3 eq 1 and M gt 2*Q;
    assert w ne 1 and w^3 eq 1 and w^2 eq R!400109;
    weights := [3^i : i in [0..d-1]];
    h := [R | 0 : i in [1..Q]];
    chi := h;
    positions := [Integers() |];
    a := K!1;
    for e in [0..Q-2] do
        pos := coordinate_index(a, weights);
        Append(~positions, pos);
        chi[pos] := e mod 2 eq 0 select R!1 else R!(-1);
        preimage := coordinate_index(-(a^3+a^2), weights);
        h[preimage] +:= 1;
        a *:= alpha;
    end for;
    h[1] +:= 1;
    assert h[1] eq 2 and &+h eq R!Q;
    radix3(~h, w);
    radix3(~chi, w);
    h := [h[i]*chi[i] : i in [1..Q]];
    delete chi;
    radix3(~h, w^2);
    // Normalize once, then recover the unique integer in [-Q,Q].
    inverse := (R!Q)^-1;
    centered := [Integers() | Integers()!(v*inverse) : v in h];
    centered := [v gt (M-1) div 2 select v-M else v : v in centered];
    assert &and[Abs(v) le Q : v in centered];
    return [-centered[pos] : pos in positions], -centered[1];
end function;

function ordinary2(K, alpha, d)
    Q := #K;
    weights := [2^i : i in [0..d-1]];
    h := [Integers() | 0 : i in [1..Q]];
    a := K!1;
    for e in [0..Q-2] do
        h[coordinate_index(a^-2, weights)] := IsZero(Trace(a)) select 1 else -1;
        a *:= alpha;
    end for;
    stride := 1;
    while stride lt Q do
        for block in [0..(Q div (2*stride))-1], j in [0..stride-1] do
            pos := block*2*stride+j+1;
            x := h[pos]; y := h[pos+stride];
            h[pos] := x+y;
            h[pos+stride] := x-y;
        end for;
        stride *:= 2;
    end while;
    ordinary := [Integers() |];
    basis := [alpha^j : j in [0..d-1]];
    a := K!1;
    for e in [0..Q-2] do
        dual := 1+&+[Integers()!Trace(a*basis[j])*weights[j] : j in [1..d]];
        Append(~ordinary, -h[dual]);
        a *:= alpha;
    end for;
    return ordinary, -h[1];
end function;

intrinsic FibrationTraceTable(p::RngIntElt, d::RngIntElt) -> Rec
{The exact ordinary trace table over GF(p^d), for p in [2,3] and 1 <= d <= 12.
The record retains Polynomial, Field, Generator, Ordinary and NodalTrace.
Ordinary[e+1] is the trace at parameter Generator^e; zero has NodalTrace}
    require p in [2, 3] and d ge 1 and d le 12:
        "p must be 2 or 3 and degree must be between 1 and 12";
    f := trace_polynomial(p, d);
    K<alpha> := ext<BaseRing(Parent(f)) | f>;
    // Magma names the degree-1 generator 1 even for the polynomial X+1 over F_3.
    if d eq 1 then
        alpha := K!(-Coefficient(f, 0));
    end if;
    assert Evaluate(f, alpha) eq 0 and Order(alpha) eq p^d-1;
    if p eq 3 then
        ordinary, nodal := ordinary3(K, alpha, d);
    else
        ordinary, nodal := ordinary2(K, alpha, d);
    end if;
    assert nodal eq 1;
    return rec<TRACE_TABLE_FORMAT | Polynomial := f, Field := K,
        Generator := alpha, Ordinary := ordinary, NodalTrace := nodal>;
end intrinsic;
