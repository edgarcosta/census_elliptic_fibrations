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
    T, mp := TorsionSubgroup(E);
    return Invariants(T), [E | mp(g) : g in T];
end intrinsic;
