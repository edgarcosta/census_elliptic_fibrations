// The action and normalization correspondence are proved in docs/mathematics.md.

import "tests/support/_weierstrass_change.m": equation_change;

// Checks the base-action order, the cross term in composition, and normalization.
procedure test_burnside_changes(p)
    k := GF(p);
    R<s, t> := PolynomialRing(k, 2);
    degrees := [2, 4, 6, 8, 12];
    A := [s^d+s^(d-1)*t+t^d : d in degrees];
    base_change := func<B, g | [Evaluate(f,
        [g[1, 1]*s+g[1, 2]*t, g[2, 1]*s+g[2, 2]*t]) : f in B]>;
    g := Matrix(k, 2, 2, [1, 1, 0, 1]);
    h := Matrix(k, 2, 2, [0, 1, 1, 0]);
    r1 := s^4+s*t^3;
    r2 := t^4;
    q1 := p eq 2 select s^2 else R!0;
    q2 := p eq 2 select t^2 else R!0;
    v1 := p eq 2 select s^6 else R!0;
    v2 := p eq 2 select t^6 else R!0;
    B := equation_change(base_change(A, g), k!1, r1, q1, v1);
    sequential := equation_change(base_change(B, h), k!1, r2, q2, v2);
    rh, qh, vh := Explode(base_change([r1, q1, v1], h));
    assert sequential eq equation_change(base_change(A, g*h), k!1,
        r2+rh, q2+qh, v2+qh*r2+vh);
    ri, qi, vi := Explode(base_change([r1, q1, v1], g^-1));
    assert A eq equation_change(base_change(B, g^-1), k!1, -ri, -qi, -vi+qi*ri);
    if p eq 2 then
        assert equation_change(A, k!1, A[2], R!0, R!0)[2] eq 0;
    else
        B := equation_change(A, k!1, R!0, -A[1]/2, -A[3]/2);
        assert B[1] eq 0 and B[3] eq 0;
        S<x, y> := PolynomialRing(R, 2);
        equation := y^2-x^3-B[2]*x^2-B[4]*x-B[5];
        assert Evaluate(equation, [x, -y]) eq equation;
        assert base_change(B, -IdentityMatrix(k, 2)) eq B;
    end if;
end procedure;

function burnside_oracle(p)
    k := GF(p);
    R<s, t> := PolynomialRing(k, 2);
    degrees := p eq 2 select [2, 4, 6, 8, 12] else [4, 8, 12];
    n := &+[d+1 : d in degrees];
    basis := [[s^(d-j)*t^j : j in [0..d]] : d in degrees];
    coefficients := func<forms | Vector(k, &cat[
        [MonomialCoefficient(forms[i], m) : m in basis[i]]
        : i in [1..#degrees]])>;
    forms := func<d | [&+[v[j+1]*s^(d-j)*t^j : j in [0..d]]
        : v in VectorSpace(k, d+1)]>;
    P := PolynomialRing(R, #degrees);
    model := p eq 2 select [P.i : i in [1..5]] else [P | 0, P.1, 0, P.2, P.3];
    indices := p eq 2 select [1..5] else [2, 4, 5];
    base_actions := [];
    for g in GL(2, k), u in k do
        if u eq 0 then continue; end if;
        coordinates := [g[1, 1]*s+g[1, 2]*t, g[2, 1]*s+g[2, 2]*t];
        rows := [];
        for i in [1..#degrees], m in basis[i] do
            image := [R | 0 : d in degrees];
            image[i] := u^(-(degrees[i] div 2))*Evaluate(m, coordinates);
            Append(~rows, coefficients(image));
        end for;
        action := Matrix(rows);
        sample := [R | &+b : b in basis];
        expected := [u^(-(degrees[i] div 2))*Evaluate(sample[i], coordinates)
            : i in [1..#degrees]];
        assert coefficients(sample)*action eq coefficients(expected);
        Append(~base_actions, action);
    end for;
    quartics := forms(4);
    quadratics := p eq 2 select forms(2) else [R!0];
    sextics := p eq 2 select forms(6) else [R!0];
    histogram := [Integers() | 0 : d in [0..n]];
    empty := 0;
    identity := IdentityMatrix(k, n);
    for r in quartics, q in quadratics, v in sextics do
        changed := equation_change(model, k!1, r, q, v)[indices];
        assert &and[TotalDegree(f) le 1 : f in changed];
        offset := coefficients([R | MonomialCoefficient(f, P!1) : f in changed]);
        rows := [];
        for j in [1..#degrees] do
            multipliers := [R | MonomialCoefficient(f, P.j) : f in changed];
            for m in basis[j] do
                Append(~rows, coefficients([m*f : f in multipliers]));
            end for;
        end for;
        translation := Matrix(rows);
        if r eq s^4 and q eq quadratics[#quadratics] and v eq sextics[#sextics] then
            sample := p eq 2 select [R | &+b : b in basis]
                else [R | 0, &+basis[1], 0, &+basis[2], &+basis[3]];
            expected := equation_change(sample, k!1, r, q, v)[indices];
            assert coefficients(sample[indices])*translation+offset eq coefficients(expected);
        end if;
        for base in base_actions do
            linear := base*translation-identity;
            if IsConsistent(linear, -offset) then
                dimension := n-Rank(linear);
                histogram[dimension+1] +:= 1;
            else
                empty +:= 1;
            end if;
        end for;
    end for;
    order := #base_actions*#quartics*#quadratics*#sextics;
    assert &+histogram+empty eq order;
    numerator := &+[histogram[d+1]*p^d : d in [0..n]];
    assert numerator mod order eq 0;
    return numerator div order, order, histogram, empty;
end function;

// Pins the whole-space orbit total against the public census entry point.
procedure test_burnside(p)
    total, order, histogram, empty := burnside_oracle(p);
    assert order eq (p eq 2 select 196608 else 23328);
    assert total eq (p eq 2 select 1402468 else 1307814237);
    expected := p eq 2 select
        [<6, 4096>, <9, 30720>, <10, 1536>, <13, 2048>, <14, 6144>,
        <15, 2880>, <18, 224>, <19, 576>, <21, 192>, <27, 1023>, <37, 1>]
        else [<3, 3888>, <5, 1728>, <7, 10368>, <9, 1728>, <10, 864>,
        <13, 968>, <15, 324>, <27, 4>];
    assert [<d-1, histogram[d]> : d in [1..#histogram] | histogram[d] ne 0] eq expected;
    assert empty eq (p eq 2 select 147168 else 3456);
    R<s, t> := PolynomialRing(GF(p), 2);
    assert total eq NumberOfFibrations(R);
end procedure;

test_burnside_changes(2);
test_burnside_changes(3);
test_burnside(2);
test_burnside(3);
