function tate_homogenize(R, a, d)
    return &+[R | Coefficient(a, j)*R.1^j*R.2^(d-j) : j in [0..d]];
end function;

function tate_homogeneous_model(R, W)
    degrees := [2, 4, 6, 8, 12];
    assert &and[Degree(W[i]) le degrees[i] : i in [1..5]];
    return [R | tate_homogenize(R, W[i], degrees[i]) : i in [1..5]];
end function;

function tate_valuation(f, pi)
    assert f ne 0;
    n := 0;
    while f mod pi eq 0 do
        f := f div pi;
        n +:= 1;
    end while;
    return n;
end function;

function tate_place_key(pi)
    R := Parent(pi);
    p := Characteristic(R);
    d := TotalDegree(pi);
    return &+[Integers()!MonomialCoefficient(pi, R.1^(d-j)*R.2^j)*p^(d-j)
        : j in [0..d]];
end function;

// LocalInformation returns the minimal discriminant valuation in info[2].
// https://magma.maths.usyd.edu.au/magma/handbook/text/1575
function tate_nonminimal(E, pi)
    F := BaseRing(E);
    U := Parent(pi);
    a := aInvariants(E);
    assert &and[Denominator(c) eq 1 : c in a];
    oldv := tate_valuation(U!Discriminant(E), pi);
    info, Emin := LocalInformation(E, F!pi);
    delta := Discriminant(Emin);
    newv := tate_valuation(U!Numerator(delta), pi)
        - tate_valuation(U!Denominator(delta), pi);
    assert newv eq info[2] and newv ge 0;
    assert oldv ge newv and (oldv-newv) mod 12 eq 0;
    return oldv gt newv, info;
end function;

function tate_bad_places(W)
    R := Universe(W);
    k := BaseRing(R);
    U<z> := PolynomialRing(k);
    F<T> := FunctionField(k);
    E := EllipticCurve([F | Evaluate(a, [T, F!1]) : a in W]);
    delta := U!Discriminant(E);
    assert delta ne 0 and Degree(delta) le 24;
    primes := [U | z+a : a in k];
    primes cat:= [U | z^2+a*z+b : a, b in k | IsIrreducible(z^2+a*z+b)];
    for factor in Factorization(delta) do
        if factor[1] notin primes then Append(~primes, factor[1]); end if;
    end for;
    bad := [R |];
    for pi in primes do
        nonminimal, info := tate_nonminimal(E, pi);
        if nonminimal then
            Append(~bad, tate_homogenize(R, pi, Degree(pi)));
        end if;
    end for;
    // The O(2i) trivialization makes every coefficient integral at infinity.
    Einf := EllipticCurve([F | Evaluate(a, [F!1, T]) : a in W]);
    assert tate_valuation(U!Discriminant(Einf), z) eq 24-Degree(delta);
    nonminimal, info := tate_nonminimal(Einf, z);
    if nonminimal then Append(~bad, R.2); end if;
    Sort(~bad, func<a, b | TotalDegree(a) ne TotalDegree(b)
        select TotalDegree(a)-TotalDegree(b) else tate_place_key(a)-tate_place_key(b)>);
    // A nonminimal place divides Delta; all other finite places have valuation 0.
    return bad;
end function;

function check_smooth_tate_certificate(W)
    bad := tate_bad_places(W);
    ok, certificate := IsEverywhereMinimal(W);
    error if certificate ne <"finite", bad> or ok ne IsEmpty(bad),
        Sprintf("Tate certificate mismatch for %o: expected %o, got %o", W, bad, certificate);
    kind, minimality := FibrationProperties(W);
    assert kind eq "smooth";
    assert minimality eq (ok select "everywhere minimal" else "nonminimal somewhere");
    return bad;
end function;

// Pins complete certificates for smooth fixtures already used by minimality tests.
procedure test_smooth_tate_fixtures(p)
    R<s, t> := PolynomialRing(GF(p), 2);
    U<z> := PolynomialRing(GF(p));
    F<T> := FunctionField(GF(p));
    if p eq 2 then
        W := [R | s*t, 0, 0, 0, s^6*t^6];
        pi := s^2+s*t+t^2;
        V := [R | pi, 0, 0, 0, pi^6];
    else
        W := [R | 0, s^2*t^2, 0, 0, s^6*t^6];
        pi := s^2+t^2;
        V := [R | 0, pi^2, 0, 0, pi^6];
    end if;
    assert check_smooth_tate_certificate(W) eq [R | t, s];
    assert check_smooth_tate_certificate(V) eq [R | pi];
    pi := s*(s+t);
    W := p eq 2 select [R | pi, 0, 0, 0, pi^6]
        else [R | 0, pi^2, 0, 0, pi^6];
    assert check_smooth_tate_certificate(W) eq [R | s, s+t];
    primes := [U | z+a : a in GF(p)];
    primes cat:= [U | z^2+a*z+b : a, b in GF(p) | IsIrreducible(z^2+a*z+b)];
    assert #primes eq (p eq 2 select 3 else 6);
    for pi in primes do
        d := Degree(pi);
        if p eq 2 then
            q := z^(d-1)+1;
            r := z^(2*d-1)+1;
            v := z^(3*d-1)+1;
            A := [U | pi, r+q*pi+q^2, r*pi, r^2+(v+r*q)*pi,
                pi^6+r^3+r*v*pi+v^2];
            wild := [U | 0, 0, pi^3, 0, pi];
        else
            r := z^(2*d-1)+1;
            A := [U | 0, pi^2, 0, -pi^2*r, pi^6+pi^2*r^2+r^3];
            wild := [U | 0, 0, 0, pi^4, pi];
        end if;
        B := tate_homogeneous_model(R, A);
        // A degree-one scaling also leaves a nonminimal point at infinity.
        expected := [R | tate_homogenize(R, pi, d)];
        if d eq 1 then expected := [R | t] cat expected; end if;
        assert check_smooth_tate_certificate(B) eq expected;
        Bswap := [R | Evaluate(a, [t, s]) : a in B];
        swapped := [R | Evaluate(a, [t, s]) : a in expected];
        swapped := [R | a/MonomialCoefficient(a, LeadingMonomial(a)) : a in swapped];
        Sort(~swapped, func<a, b | TotalDegree(a) ne TotalDegree(b)
            select TotalDegree(a)-TotalDegree(b) else tate_place_key(a)-tate_place_key(b)>);
        assert check_smooth_tate_certificate(Bswap) eq swapped;

        E := EllipticCurve([F | a : a in wild]);
        nonminimal, info := tate_nonminimal(E, pi);
        // Additive reduction has tame conductor exponent 2; the excess is wild.
        assert not nonminimal and info[2] eq 12 and info[3] gt 2;
        B := tate_homogeneous_model(R, wild);
        expected := d eq 1 select [R | t] else [R |];
        assert check_smooth_tate_certificate(B) eq expected;
    end for;
end procedure;

// A missed place or a false positive must disagree with local Tate valuations.
procedure test_smooth_tate_sample(p)
    R<s, t> := PolynomialRing(GF(p), 2);
    F<T> := FunctionField(GF(p));
    n := p eq 2 select 21002 else 31003;
    bound := p eq 2 select 2^32 else 3^27;
    smooth := 0;
    minimal := 0;
    while smooth lt 256 do
        n := (1664525*n+1013904223) mod bound;
        W := aInvariantsFromKey(R, n);
        // Tate's algorithm applies only when the generic fibre is smooth.
        a := [F | Evaluate(c, [T, F!1]) : c in W];
        b2 := a[1]^2+4*a[2];
        b4 := a[1]*a[3]+2*a[4];
        b6 := a[3]^2+4*a[5];
        b8 := a[1]^2*a[5]+4*a[2]*a[5]-a[1]*a[3]*a[4]+a[2]*a[3]^2-a[4]^2;
        delta := -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6;
        if delta eq 0 then continue; end if;
        bad := check_smooth_tate_certificate(W);
        smooth +:= 1;
        if IsEmpty(bad) then minimal +:= 1; end if;
    end while;
    assert minimal gt 0;
end procedure;
