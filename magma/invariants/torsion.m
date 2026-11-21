intrinsic FibrationTorsion(E::CrvEll[FldFunRat]) -> SeqEnum, SeqEnum
{The torsion invariant factors of E and all its torsion points, with the identity
first, for an everywhere minimal elliptic K3 polynomial equation over GF(2)(t)
or GF(3)(t). See also FibrationLocalData}
    F := BaseRing(E);
    k := BaseRing(F);
    require Type(k) eq FldFin and #k in [2, 3]:
        "constant field must be GF(2) or GF(3)";
    R := PolynomialRing(k, 2);
    W := aInvariants(E, R);
    kind, minimality := FibrationProperties(W);
    require kind eq "smooth" and minimality eq "everywhere minimal":
        "equation must define an everywhere minimal elliptic K3";
    T, mp := TorsionSubgroup(E);
    return Invariants(T), [E | mp(g) : g in T];
end intrinsic;
