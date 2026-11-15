INVARIANTS_FORMAT := recformat<Rank, GeometricRank, PicardNumber,
    GeometricPicardNumber, P2, ZetaSign, PointCounts, BrauerHeight,
    ArtinTateProduct, NeronSeveriSplittingDegree, JDegree,
    JInseparabilityExponent, IsIsotrivial, ShaRegulatorProduct>;

function quotient_at_one(L, q)
    P := Parent(L);
    factor := 1-q*P.1;
    rank := 0;
    while Evaluate(L, 1/q) eq 0 do
        L := ExactQuotient(L, factor);
        rank +:= 1;
    end while;
    return L, rank;
end function;

intrinsic FibrationInvariants(E::CrvEll[FldFunRat], L::RngUPolElt,
    fibres::SeqEnum, torsion::SeqEnum) -> Rec
{The exact derived K3 invariants of E from its complete L-polynomial L,
FibrationLocalData fibre records and torsion invariant factors over a finite
constant field. The supplied data are assumed correct for E.
See also FibrationLPolynomial, FibrationTorsion}
    k := BaseRing(BaseRing(E));
    require Type(k) eq FldFin: "constant field must be finite";
    q := #k;
    p := Characteristic(k);
    f := Degree(k);
    P<T> := PolynomialRing(Rationals());
    require &and[IsCoercible(Integers(), c) : c in Coefficients(L)]
        and Coefficient(L, 0) eq 1:
        "L must have integer coefficients and constant coefficient 1";
    require &and[Type(n) eq RngIntElt and n gt 1 : n in torsion]
        and &and[torsion[i+1] mod torsion[i] eq 0 : i in [1..#torsion-1]]:
        "torsion must contain invariant factors greater than 1 in divisibility order";
    L := P!L;
    cycles := &cat[[Integers() | n : n in r`ComponentCycles] : r in fibres];
    P2 := L*(1-q*T)^2;
    for n in cycles do P2 *:= 1-(q*T)^n; end for;
    require Degree(P2) eq 22: "L and component cycles must give P2 of degree 22";
    Q, rank := quotient_at_one(L, q);
    G, rho := quotient_at_one(P2, q);
    assert rho eq 2+rank+#cycles;
    geometric_rank := 0;
    splitting := LCM([1] cat cycles);
    for factor in Factorization(Evaluate(L, T/q)) do
        cyclotomic, order := IsCyclotomicPolynomial(factor[1]);
        if cyclotomic then
            geometric_rank +:= Degree(factor[1])*factor[2];
            splitting := LCM(splitting, order);
        end if;
    end for;
    slope := Min([Rationals() | 1] cat [Valuation(Integers()!Coefficient(L, i), p)/(f*i)
        : i in [1..Degree(L)] | Coefficient(L, i) ne 0]);
    if slope eq 1 then
        height := Infinity();
    else
        height := Integers()!(1/(1-slope));
        assert height in [1..10];
    end if;
    coefficients := Coefficients(P2);
    sums := [Integers() |];
    for n in [1..5] do
        Append(~sums, Integers()!(-n*coefficients[n+1]
            - &+[Rationals() | coefficients[j+1]*sums[n-j] : j in [1..n-1]]));
    end for;
    artin_tate := q*Evaluate(G, 1/q);
    assert IsIntegral(artin_tate) and artin_tate gt 0;
    tamagawa := &*[Integers() | r`Tamagawa : r in fibres];
    torsion_order := &*[Integers() | n : n in torsion];
    sha_regulator := q*Evaluate(Q, 1/q)*torsion_order^2/tamagawa;
    j := jInvariant(E);
    j_degree := Max(Degree(Numerator(j)), Degree(Denominator(j)));
    inseparability := 0;
    if j_degree gt 0 then
        while Derivative(j) eq 0 do
            j := Root(j, p);
            inseparability +:= 1;
        end while;
    end if;
    return rec<INVARIANTS_FORMAT | Rank := rank, GeometricRank := geometric_rank,
        PicardNumber := rho, GeometricPicardNumber := 2+geometric_rank
            + &+([0] cat cycles), P2 := PolynomialRing(Integers())!P2,
        ZetaSign := Integers()!(LeadingCoefficient(P2)/q^22),
        PointCounts := [1+q^(2*n)+sums[n] : n in [1..5]], BrauerHeight := height,
        ArtinTateProduct := Integers()!artin_tate,
        NeronSeveriSplittingDegree := splitting, JDegree := j_degree,
        JInseparabilityExponent := inseparability, IsIsotrivial := j_degree eq 0,
        ShaRegulatorProduct := sha_regulator>;
end intrinsic;

intrinsic FibrationInvariants(E::CrvEll[FldFunRat]) -> Rec
{The exact derived invariants of an everywhere minimal elliptic K3 polynomial
model E over GF(2)(t) or GF(3)(t). See also FibrationLPolynomial,
FibrationLocalData, FibrationTorsion}
    L := FibrationLPolynomial(E);
    fibres := FibrationLocalData(E);
    torsion := FibrationTorsion(E);
    return FibrationInvariants(E, L`Polynomial, fibres, torsion);
end intrinsic;
