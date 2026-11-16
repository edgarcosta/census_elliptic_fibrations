SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/_native_derived_fields.m": derived_fields;

procedure write_constant_field_controls(path)
    file := Open(path, "w");
    delete file;
    F<t> := FunctionField(GF(2^12));
    E := EllipticCurve([F | 0, 0, 1, t^6, 0]);
    P<T> := PolynomialRing(Integers());
    R := recformat<ComponentCycles, Tamagawa>;
    rows := [rec<R | ComponentCycles := [1 : i in [1..13]], Tamagawa := 4>];
    result := FibrationInvariants(E, (1-4096*T)^7, rows, []);
    assert LeadingCoefficient(result`P2) eq 2^264;
    assert result`PointCounts[5] eq 1329227995784915898268080161630978049;
    fprintf path, "%o\n", Join(derived_fields(result), ",");
    F<t> := FunctionField(GF(2));
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    P<T> := PolynomialRing(Rationals());
    M := CompanionMatrix(P!Reverse(Coefficients(LFunction(E))));
    extended := Determinant(ScalarMatrix(Nrows(M), P!1)-T*ChangeRing(M^2, P));
    F<t> := FunctionField(GF(4));
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    rows := [rec<R | ComponentCycles := [1], Tamagawa := 2>];
    result := FibrationInvariants(E, extended, rows, []);
    assert result`BrauerHeight eq 3;
    fprintf path, "%o\n", Join(derived_fields(result), ",");
end procedure;

write_constant_field_controls(referencepath);
quit;
