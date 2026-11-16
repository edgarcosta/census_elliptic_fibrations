import "tests/_weierstrass_change.m": equation_change;

function invariance_form(R, d, n)
    f := R!0;
    for j in [0..d] do
        n := (1664525*n+1013904223) mod 2^32;
        f +:= BaseRing(R)!((n div 2^16) mod Characteristic(R))*R.1^(d-j)*R.2^j;
    end for;
    return f, n;
end function;

function certificate_set(certificate)
    return certificate[1] eq "all" select certificate
        else <"finite", Seqset(certificate[2])>;
end function;

function transported_certificate(certificate, coordinates)
    if certificate[1] eq "all" then return certificate; end if;
    places := {Evaluate(pi, coordinates) : pi in certificate[2]};
    return <"finite", {pi/MonomialCoefficient(pi, LeadingMonomial(pi)) : pi in places}>;
end function;

function public_signature(W)
    kind, state := FibrationProperties(W);
    minimal, certificate := IsEverywhereMinimal(W);
    return <kind, state, minimal, certificate_set(certificate)>;
end function;

procedure check_invariance_model(W, expected)
    actual := public_signature(W);
    error if actual ne expected,
        Sprintf("Coordinate-change mismatch for %o: expected %o, got %o", W, expected, actual);
    R := Universe(W);
    key := FibrationKey(W);
    decoded := aInvariantsFromKey(R, key);
    // Normalization is a coordinate change; keys need not agree on other equations.
    normalized := Characteristic(R) eq 2
        select equation_change(W, BaseRing(R)!1, W[2], R!0, R!0)
        else equation_change(W, BaseRing(R)!1, R!0, -W[1]/2, -W[3]/2);
    error if decoded ne normalized,
        Sprintf("Decoded normalization mismatch for %o", W);
    error if public_signature(decoded) ne expected,
        Sprintf("Decoded properties or certificate changed for %o", W);
    assert FibrationKey(decoded) eq key;
end procedure;

function invariance_fixtures(R)
    s := R.1;
    t := R.2;
    if Characteristic(R) eq 2 then
        smooth := [R | s*t, 0, s^6+t^6, s^7*t, s^11*t];
        multiple := [R | s*t, 0, 0, 0, s^6*t^6];
        pi := s^2+s*t+t^2;
        quadratic := [R | pi, 0, 0, 0, pi^6];
        nodal := [R | s*t, s^4+t^4, 0, 0, 0];
        nonregular := [R | 0, 0, 0, s^6*t^2, 0];
        everywhere := [R | 0, 0, 0, s^8, t^12];
    else
        smooth := [R | 0, s*t*(s^2+t^2), 0, t^8, s^12];
        multiple := [R | 0, s^2*t^2, 0, 0, s^6*t^6];
        pi := s^2+t^2;
        quadratic := [R | 0, pi^2, 0, 0, pi^6];
        nodal := [R | 0, s*t*(s^2+t^2), 0, 0, 0];
        nonregular := [R | 0, s^4, 0, 0, 0];
        everywhere := [R | 0, 0, 0, 0, s^12];
    end if;
    return [* <smooth, "smooth", <"finite", [R |]>>,
        <multiple, "smooth", <"finite", [R | t, s]>>,
        <quadratic, "smooth", <"finite", [R | pi]>>,
        <[R | 0, 0, 0, 0, s^11*t+s*t^11], "regular cusp", <"finite", [R |]>>,
        <[R | 0, 0, 0, 0, s*t^11], "regular cusp", <"finite", [R | t]>>,
        <nodal, "nonregular", <"finite", [R |]>>,
        <nonregular, "nonregular", <"finite", [R | s]>>,
        <[R | 0, 0, 0, 0, 0], "nonregular", <"all">>,
        <everywhere, "nonregular", <"all">> *];
end function;

// Pins classification, place transport and decoded normalization in all six cells.
procedure test_coordinate_invariance(p)
    k := GF(p);
    R<s, t> := PolynomialRing(k, 2);
    n := p eq 2 select 42002 else 43003;
    fixtures := invariance_fixtures(R);
    for j in [1..8] do
        W := [R |];
        for d in [2, 4, 6, 8, 12] do
            f, n := invariance_form(R, d, n);
            Append(~W, f);
        end for;
        kind, state := FibrationProperties(W);
        minimal, certificate := IsEverywhereMinimal(W);
        Append(~fixtures, <W, kind, certificate>);
    end for;
    translations := [<R!0, R!0, R!0>];
    for j in [1..3] do
        r, n := invariance_form(R, 4, n);
        q, n := invariance_form(R, 2, n);
        v, n := invariance_form(R, 6, n);
        Append(~translations, <r, q, v>);
    end for;
    swap := Matrix(k, 2, 2, [0, 1, 1, 0]);
    shear := Matrix(k, 2, 2, [1, 1, 0, 1]);
    bases := [IdentityMatrix(k, 2), swap, shear, shear*swap];
    if p eq 3 then Append(~bases, DiagonalMatrix(k, [2, 1])); end if;
    cells := {};
    moved := false;
    rescaled := false;
    for fixture in fixtures do
        W, kind, certificate := Explode(fixture);
        minimal := certificate[1] eq "finite" and IsEmpty(certificate[2]);
        state := minimal select "everywhere minimal" else "nonminimal somewhere";
        Include(~cells, <kind, minimal>);
        check_invariance_model(W, <kind, state, minimal, certificate_set(certificate)>);
        for g in bases do
            coordinates := [g[1, 1]*s+g[1, 2]*t, g[2, 1]*s+g[2, 2]*t];
            B := [Evaluate(a, coordinates) : a in W];
            expected := transported_certificate(certificate, coordinates);
            moved or:= expected ne certificate_set(certificate);
            if certificate[1] eq "finite" then
                rescaled or:= &or[LeadingCoefficient(Evaluate(pi, coordinates)) ne 1
                    : pi in certificate[2]];
            end if;
            for translation in translations, u in k do
                if u eq 0 then continue; end if;
                r, q, v := Explode(translation);
                changed := equation_change(B, u, r, q, v);
                check_invariance_model(changed, <kind, state, minimal, expected>);
                assert equation_change(changed, u^-1, -r/u^2, -q/u,
                    (q*r-v)/u^3) eq B;
            end for;
        end for;
    end for;
    assert cells eq {<kind, minimal> : kind in ["smooth", "regular cusp", "nonregular"],
        minimal in [true, false]};
    assert moved;
    if p eq 3 then assert rescaled; end if;
end procedure;

test_coordinate_invariance(2);
test_coordinate_invariance(3);
