function geometric_correction(row)
    case row`Type:
        when "I": return Floor(row`Index^2/4)/row`Index;
        when "I*": return (row`Index+4)/4;
        when "III": return 1/2;
        when "IV": return 2/3;
        when "IV*": return 4/3;
        when "III*": return 3/2;
        else return 0;
    end case;
end function;

function zero_height_sections(E, bound)
    F := BaseRing(E);
    k := BaseRing(F);
    R<t> := PolynomialRing(k);
    S<y> := PolynomialRing(F);
    a1, a2, a3, a4, a6 := Explode([R!a : a in aInvariants(E)]);
    points := [E | E!0];
    if bound lt 0 then return points; end if;
    denominators := [R | 1] cat [t+c : c in k | bound eq 1];
    for g in denominators do
        g2 := g^2;
        g3 := g^3;
        for v in CartesianPower(k, 5+2*bound) do
            A := R![v[i] : i in [1..5+2*bound]];
            if GCD(A, g) ne 1 then continue; end if;
            rhs := A^3+a2*A^2*g2+a4*A*g2^2+a6*g3^2;
            linear := a1*A*g+a3*g3;
            for root in Roots(y^2+(F!linear)*y-(F!rhs)) do
                B := root[1];
                if Denominator(B) ne 1 or Degree(Numerator(B)) gt 6+3*bound then
                    continue;
                end if;
                infinity := Max(0, Ceiling((Degree(A)-2*Degree(g)-4)/2));
                if Degree(g)+infinity gt bound then continue; end if;
                P := E![(F!A)/g2, B/g3];
                if Height(P) eq 0 then Append(~points, P); end if;
            end for;
        end for;
    end for;
    return points;
end function;

intrinsic FibrationTorsion(E::CrvEll[FldFunRat]) -> SeqEnum, SeqEnum
{The torsion invariant factors of E and all its torsion points, with the identity
first, for an everywhere minimal elliptic K3 polynomial equation over GF(2)(t)
or GF(3)(t). See also FibrationLocalData}
    F := BaseRing(E);
    k := BaseRing(F);
    require Type(k) eq FldFin and #k in [2, 3]:
        "constant field must be GF(2) or GF(3)";
    a := aInvariants(E);
    weights := [2, 4, 6, 8, 12];
    require &and[Denominator(a[i]) eq 1 and Degree(Numerator(a[i])) le weights[i]
        : i in [1..5]]:
        "coefficients must be polynomials of degrees at most 2,4,6,8,12";
    R<s, u> := PolynomialRing(k, 2);
    W := [R | &+[R | Coefficient(Numerator(a[i]), j)*s^j*u^(weights[i]-j)
        : j in [0..Degree(Numerator(a[i]))]] : i in [1..5]];
    kind, minimality := FibrationProperties(W);
    require kind eq "smooth" and minimality eq "everywhere minimal":
        "equation must define an everywhere minimal elliptic K3";
    rows := FibrationLocalData(E);
    correction := &+[Rationals() | row`Degree*geometric_correction(row) : row in rows];
    assert correction le 6;
    // The degree-weighted correction bound includes wild fibres and infinity.
    bound := Floor((correction-4)/2);
    points := zero_height_sections(E, bound);
    orders := [Integers() |];
    for P in points do
        Q := E!0;
        for n in [1..#points] do
            Q +:= P;
            if IsId(Q) then Append(~orders, n); break; end if;
        end for;
    end for;
    assert #orders eq #points;
    // Elliptic torsion has at most two invariant factors, including p-torsion.
    exponent := LCM(orders);
    factors := [#points div exponent, exponent];
    assert factors[2] mod factors[1] eq 0;
    return [n : n in factors | n gt 1], points;
end intrinsic;
