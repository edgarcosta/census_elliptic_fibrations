function gf4096_fixture()
    F<t> := FunctionField(GF(2^12));
    E := EllipticCurve([F | 0, 0, 1, t^6, 0]);
    P<T> := PolynomialRing(Integers());
    R := recformat<ComponentCycles, Tamagawa>;
    rows := [rec<R | ComponentCycles := [1 : i in [1..13]], Tamagawa := 4>];
    return E, (1-4096*T)^7, rows;
end function;

function gf4_fixture()
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
    return E, extended, rows, L;
end function;
