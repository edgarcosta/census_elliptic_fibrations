import "../_native_torsion_samples.m": torsion_samples;

function code(f, p)
    return &+[Integers() | (Integers()!Coefficient(f, i))*p^i
        : i in [0..Degree(f)]];
end function;

function coordinates(P, p)
    if IsId(P) then return "O"; end if;
    return Join([Sprint(code(f, p)) : f in [Numerator(P[1]), Denominator(P[1]),
        Numerator(P[2]), Denominator(P[2])]], ":");
end function;

procedure write_torsion_reference(path)
    empty := Open(path, "w");
    delete empty;
    fprintf path, "q,key,torsion,torsion_points\n";
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        for a in torsion_samples(F) do
            W := aInvariantsFromKey(R, FibrationKey(aInvariants(EllipticCurve(a), R)));
            for chart in [W, [R | Evaluate(f, [u, s]) : f in W]] do
                E := EllipticCurve(chart, F);
                factors, points := FibrationTorsion(E);
                labels := [coordinates(P, p) : P in points];
                fprintf path, "%o,%o,%o,%o\n", p, FibrationKey(chart),
                    IsEmpty(factors) select "[]" else Join([Sprint(n) : n in factors], " "),
                    Join(labels, ";");
            end for;
        end for;
    end for;
end procedure;
