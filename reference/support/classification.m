function residue_polynomials(U, n)
    p := #BaseRing(U);
    z := U.1;
    return [ U | &+[ U | ((a div p^j) mod p)*z^j : j in [0..n-1] ] : a in [0..p^n-1] ];
end function;

function reference_local_nonminimal(W, pi)
    a1,a2,a3,a4,a6 := Explode(W);
    U := Parent(pi);
    if Characteristic(U) eq 2 then
        if a1 mod pi ne 0 then return false; end if;
        for q in residue_polynomials(U,Degree(pi)) do
            r := (a2+q*a1+q^2) mod pi^2;
            if (a3+r*a1) mod pi^3 ne 0 then continue; end if;
            for v in residue_polynomials(U,3*Degree(pi)) do
                b4 := a4+q*a3+(v+r*q)*a1+r^2;
                b6 := a6+r*a4+r^2*a2+r^3+v*a3+r*v*a1+v^2;
                if b4 mod pi^4 eq 0 and b6 mod pi^6 eq 0 then return true; end if;
            end for;
        end for;
    else
        C := a2+a1^2/4;
        A := a4+a1*a3/2;
        B := a6+a3^2/4;
        if C mod pi^2 ne 0 then return false; end if;
        for r in residue_polynomials(U,2*Degree(pi)) do
            if (A-C*r) mod pi^4 eq 0 and (B+A*r+C*r^2+r^3) mod pi^6 eq 0 then return true; end if;
        end for;
    end if;
    return false;
end function;

function factored_power(f, n)
    return IsZero(f) or &and[ x[2] mod n eq 0 : x in Factorization(f) ];
end function;

function reference_properties(W)
    R := Universe(W);
    p := Characteristic(R);
    U<z> := PolynomialRing(BaseRing(R));
    a1,a2,a3,a4,a6 := Explode([ U | Evaluate(a,[z,U!1]) : a in W ]);
    b2 := a1^2+4*a2;
    b4 := a1*a3+2*a4;
    b6 := a3^2+4*a6;
    b8 := a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
    delta := -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6;
    if delta ne 0 then
        kind := "smooth";
    elif p eq 2 then
        kind := a1 eq 0 and a3 eq 0 and
            (not factored_power(a4+a2^2,2) or not factored_power(a6+a2*a4,2))
            select "regular cusp" else "nonregular";
    else
        kind := a2+a1^2/4 eq 0 and a4+a1*a3/2 eq 0 and not factored_power(a6+a3^2/4,3)
            select "regular cusp" else "nonregular";
    end if;
    allplaces := p eq 2 select a1 eq 0 and a3 eq 0
        and factored_power(a4+a2^2,4) and factored_power(a6+a2*a4,2)
        else a2+a1^2/4 eq 0 and a4+a1*a3/2 eq 0 and factored_power(a6+a3^2/4,3);
    // Both charts are needed: a model can be nonminimal only at infinity.
    // Quadratic places are tested on cusps too, independently of the place bound.
    primes := [ U | ];
    for d in [1,2] do
        for f in residue_polynomials(U,d) do
            pi := z^d+f;
            if IsIrreducible(pi) then Append(~primes,pi); end if;
        end for;
    end for;
    bad := [ R | ];
    affine := [ U | Evaluate(a,[z,U!1]) : a in W ];
    for pi in primes do
        if reference_local_nonminimal(affine,pi) then
            Append(~bad,&+[ R | Coefficient(pi,j)*R.1^j*R.2^(Degree(pi)-j) : j in [0..Degree(pi)] ]);
        end if;
    end for;
    if reference_local_nonminimal([ U | Evaluate(a,[U!1,z]) : a in W ],z) then Append(~bad,R.2); end if;
    if allplaces then
        assert #bad eq #primes+1;
        cert := <"all">;
    else
        key := func< f | &+[ Integers()!MonomialCoefficient(f,R.1^(TotalDegree(f)-j)*R.2^j)*p^(TotalDegree(f)-j) : j in [0..TotalDegree(f)] ] >;
        Sort(~bad,func< f,g | TotalDegree(f) eq TotalDegree(g) select key(f)-key(g) else TotalDegree(f)-TotalDegree(g) >);
        cert := <"finite",bad>;
    end if;
    return kind,#bad eq 0 select "everywhere minimal" else "nonminimal somewhere",cert;
end function;

procedure check_reference_properties(W)
    expected_kind,expected_min,expected_cert := reference_properties(W);
    kind,minimum := FibrationProperties(W);
    ok,certificate := IsEverywhereMinimal(W);
    error if kind ne expected_kind or minimum ne expected_min or certificate ne expected_cert
        or ok ne (minimum eq "everywhere minimal"),Sprintf("Reference disagreement for %o: %o versus %o",W,<kind,minimum,certificate>,<expected_kind,expected_min,expected_cert>);
end procedure;
