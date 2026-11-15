SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/_support.m": expect_error;

// Independent Magma algorithms pin ranks, cohomology, heights and both products.
procedure check_invariants(E, rank, height, splitting, isotrivial)
    p := #BaseRing(BaseRing(E));
    result := FibrationInvariants(E);
    P<T> := PolynomialRing(Rationals());
    L := P!LFunction(E);
    information := AnalyticInformation(E, L);
    M := FrobeniusActionOnTrivialLattice(E);
    P2 := L*Determinant(ScalarMatrix(Nrows(M), P!1)-p*T*ChangeRing(M, P));
    rho, product := ArtinTateFormula(P!Reverse(Coefficients(P2)), p, 1);
    assert result`Rank eq rank and rank eq information[1];
    assert result`GeometricRank eq information[2];
    assert result`P2 eq P2 and Degree(P2) eq 22;
    assert result`PicardNumber eq rho;
    assert result`GeometricPicardNumber eq 2+information[2]
        + &+[Integers() | Degree(li[1])*(NumberOfComponents(li[5])-1)
            : li in LocalInformation(E)];
    assert result`ShaRegulatorProduct eq information[3];
    assert result`ArtinTateProduct eq product;
    slope := Slopes(NewtonPolygon(ChangeRing(P2, pAdicRing(p, 100))))[1];
    assert result`BrauerHeight eq height;
    assert height eq (slope eq 1 select Infinity() else 1/(1-slope));
    assert result`NeronSeveriSplittingDegree eq splitting;
    assert result`ZetaSign eq LeadingCoefficient(P2)/p^22;
    S<z> := PowerSeriesRing(Rationals(), 7);
    Z := 1/((1-z)*(S!P2)*(1-p^2*z));
    assert result`PointCounts eq [Integers()!Coefficient(z*Derivative(Z)/Z, i)
        : i in [1..5]];
    j := jInvariant(E);
    degree := Max(Degree(Numerator(j)), Degree(Denominator(j)));
    assert result`JDegree eq degree;
    assert result`IsIsotrivial eq isotrivial and isotrivial eq (degree eq 0);
    if degree gt 0 then
        assert (Derivative(j) eq 0) eq (result`JInseparabilityExponent gt 0);
        j := Root(j, p^result`JInseparabilityExponent);
        assert Derivative(j) ne 0;
    else
        assert result`JInseparabilityExponent eq 0;
    end if;
    rows := FibrationLocalData(E);
    assert rho eq 2+rank+&+[Integers() | #r`ComponentCycles : r in rows];
    assert result`GeometricPicardNumber eq 2+result`GeometricRank
        + &+[Integers() | &+([0] cat r`ComponentCycles) : r in rows];
    // A closed place contributes its Tamagawa number once, even in degree > 1.
    Q := L;
    for i in [1..rank] do Q := ExactQuotient(Q, 1-p*T); end for;
    torsion := TorsionSubgroup(E);
    C := &*[Integers() | r`Tamagawa : r in rows];
    assert result`ShaRegulatorProduct eq p*Evaluate(Q, 1/p)*#torsion^2/C;
    f := P!1;
    g := P!1;
    for r in rows do
        for n in r`ComponentCycles do
            f *:= 1-(p*T)^n;
            g *:= 1-p*T;
        end for;
    end for;
    ratio := Evaluate(ExactQuotient(f, g), 1/p);
    assert result`ArtinTateProduct eq result`ShaRegulatorProduct*C*ratio/#torsion^2;
end procedure;

// Rank zero can have nontrivial torsion; supersingular height is infinity.
procedure test_fixed_derived_invariants()
    F<t> := FunctionField(GF(2));
    check_invariants(EllipticCurve([F | 0, 0, 1, t^6, 0]), 2, Infinity(), 6, true);
    check_invariants(EllipticCurve([F | 0, 0, 1, t^7, 0]), 1, 3, 12, true);
    check_invariants(EllipticCurve([F | t, 0, t^6+1, 0, 0]), 1, 1, 2, false);
    check_invariants(EllipticCurve([F | t^2, 0, 1, 0, 0]), 0, Infinity(), 2, false);
    F<t> := FunctionField(GF(3));
    check_invariants(EllipticCurve([F | 0, 0, 0, t^8, t+1]), 2, 8, 2, true);
    check_invariants(EllipticCurve([F | 0, t^3+t, 0, t^8+1, 0]), 2, 4, 4, false);
    check_invariants(EllipticCurve([F | 0, 0, 0, 1, t^7]), 2, Infinity(), 6, true);
    check_invariants(EllipticCurve([F | 0, t^4+1, 0, t^4, 0]),
        1, Infinity(), 2, false);
end procedure;

// Constant-field extension preserves finite formal height and normalizes v_q(q)=1.
procedure test_height_over_a_nonprime_field()
    F<t> := FunctionField(GF(2));
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    P<T> := PolynomialRing(Rationals());
    L := P!LFunction(E);
    M := CompanionMatrix(P!Reverse(Coefficients(L)));
    extended := Determinant(ScalarMatrix(Nrows(M), P!1)-T*ChangeRing(M^2, P));
    F<t> := FunctionField(GF(4));
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    R := recformat<ComponentCycles, Tamagawa>;
    rows := [rec<R | ComponentCycles := [1], Tamagawa := 2>];
    result := FibrationInvariants(E, extended, rows, []);
    assert result`BrauerHeight eq 3;
    P2 := L*(1-2*T)^3;
    M := CompanionMatrix(P!Reverse(Coefficients(P2)));
    assert result`PointCounts eq [1+4^(2*n)+Trace(M^(2*n)) : n in [1..5]];
end procedure;

// Over GF(2^12), all 22 eigenvalues of this supersingular K3 equal q.
procedure test_large_exact_integers()
    q := 2^12;
    F<t> := FunctionField(GF(q));
    E := EllipticCurve([F | 0, 0, 1, t^6, 0]);
    P<T> := PolynomialRing(Integers());
    R := recformat<ComponentCycles, Tamagawa>;
    rows := [rec<R | ComponentCycles := [1 : i in [1..13]], Tamagawa := 4>];
    result := FibrationInvariants(E, (1-q*T)^7, rows, []);
    assert result`P2 eq (1-q*T)^22;
    assert result`PointCounts eq [1+q^(2*n)+22*q^n : n in [1..5]];
    assert result`PointCounts[5] gt 2^53 and IsOdd(result`PointCounts[5]);
    assert LeadingCoefficient(result`P2) eq 2^264;
    assert result`Rank eq 7 and result`GeometricRank eq 7;
    assert result`PicardNumber eq 22 and result`GeometricPicardNumber eq 22;
    assert result`BrauerHeight eq Infinity();
    assert result`ArtinTateProduct eq q and result`ShaRegulatorProduct eq q/4;
    assert result`NeronSeveriSplittingDegree eq 1 and result`ZetaSign eq 1;
end procedure;

// Supplied polynomials, torsion, and component data must have the stated shape.
procedure test_supplied_data_domain()
    F<t> := FunctionField(GF(2));
    E := EllipticCurve([F | t^2, 0, 1, 0, 0]);
    rows := FibrationLocalData(E);
    P<T> := PolynomialRing(Rationals());
    for L in [P | 2, 1+T/2] do
        expect_error(procedure() result := FibrationInvariants(E, L, rows, [6]);
            end procedure, "L must have integer coefficients and constant coefficient 1");
    end for;
    expect_error(procedure() result := FibrationInvariants(E, P!1, rows, [4, 2]);
        end procedure, "torsion must contain invariant factors");
    expect_error(procedure() result := FibrationInvariants(E, 1+T, rows, [6]);
        end procedure, "L and component cycles must give P2 of degree 22");
end procedure;

test_fixed_derived_invariants();
test_height_over_a_nonprime_field();
test_large_exact_integers();
test_supplied_data_domain();
Write(resultpath, "DERIVED_MAGMA_OK");
quit;
