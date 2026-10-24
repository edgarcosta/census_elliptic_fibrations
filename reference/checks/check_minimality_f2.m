AttachSpec("spec");
load "reference/f2/elliptic_fibration_orbits_f2.m";
load "reference/support/classification.m";

function polynomial_valuation(f,pi)
    assert f ne 0;
    n := 0;
    while f mod pi eq 0 do
        f := f div pi;
        n +:= 1;
    end while;
    return n;
end function;

function homogeneous_model(R,W)
    degrees := [2,4,6,8,12];
    assert &and[ Degree(W[i]) le degrees[i] : i in [1..5] ];
    return [ R | &+[ R | Coefficient(W[i],j)*R.1^j*R.2^(degrees[i]-j)
        : j in [0..degrees[i]] ] : i in [1..5] ];
end function;

function tate_bad_places(W)
    R := Universe(W);
    U<z> := PolynomialRing(GF(2));
    F<T> := FunctionField(GF(2));
    affine := [ U | Evaluate(a,[z,U!1]) : a in W ];
    E := EllipticCurve([ F | a : a in affine ]);
    delta := U!Discriminant(E);
    bad := [ R | ];
    for pi in [z,z+1,z^2+z+1] do
        info,Emin := LocalInformation(E,F!pi);
        newdelta := Discriminant(Emin);
        newv := polynomial_valuation(U!Numerator(newdelta),pi)
            - polynomial_valuation(U!Denominator(newdelta),pi);
        oldv := polynomial_valuation(delta,pi);
        assert newv eq info[2] and oldv ge newv and (oldv-newv) mod 12 eq 0;
        if oldv gt newv then
            Append(~bad,&+[ R | Coefficient(pi,j)*R.1^j*R.2^(Degree(pi)-j)
                : j in [0..Degree(pi)] ]);
        end if;
    end for;
    // The equation has coefficients in O(2i), so infinity uses degree 24.
    info,Emin := LocalInformation(E,1/T);
    newdelta := Discriminant(Emin);
    newv := Degree(Denominator(newdelta))-Degree(Numerator(newdelta));
    oldv := 24-Degree(delta);
    assert newv eq info[2] and oldv ge newv and (oldv-newv) mod 12 eq 0;
    if oldv gt newv then Append(~bad,R.2); end if;
    return SequenceToSet(bad);
end function;

// Tate and direct congruences check both charts, including wild reduction.
procedure test_f2_minimality_against_tate()
    R<s,t> := PolynomialRing(GF(2),2);
    U<z> := PolynomialRing(GF(2));
    models := [];
    for pi in [z,z+1,z^2+z+1] do
        d := Degree(pi);
        q := z^(d-1)+1; r := z^(2*d-1)+1; v := z^(3*d-1)+1;
        W := [ U | pi,0,0,0,pi^6 ];
        translated := [ U | pi,r+q*pi+q^2,r*pi,r^2+(v+r*q)*pi,
            pi^6+r^3+r*v*pi+v^2 ];
        for a in [W,translated,[U | 0,0,pi^3,0,pi]] do
            B := homogeneous_model(R,a);
            Append(~models,B);
            Append(~models,[ R | Evaluate(b,[t,s]) : b in B ]);
        end for;
    end for;
    // An integer recurrence fixes the sample without changing Magma's seed.
    n := 21002;
    for i in [1..64] do
        n := (1664525*n+1013904223) mod 2^32;
        W := FibrationFromKeyF2(R,n);
        if FibrationDiscriminantF2(W) ne 0 then Append(~models,W); end if;
    end for;
    assert #models ge 64;
    minimal := 0; nonminimal := 0;
    for W in models do
        bad := tate_bad_places(W);
        check_reference_properties(W);
        ok,certificate := IsEverywhereMinimal(W);
        assert certificate[1] eq "finite" and SequenceToSet(certificate[2]) eq bad;
        assert ok eq IsEmpty(bad);
        if ok then minimal +:= 1; else nonminimal +:= 1; end if;
    end for;
    assert minimal gt 0 and nonminimal gt 0;
    printf "F2 Tate sample: smooth minimal=%o nonminimal=%o; models=%o; local checks=%o\n",
        minimal,nonminimal,#models,4*#models;
end procedure;

test_f2_minimality_against_tate();
Write(resultpath,"REFERENCE_CHECK_OK");
quit 0;
