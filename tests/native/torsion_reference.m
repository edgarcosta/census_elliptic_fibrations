SetQuitOnError(true);
AttachSpec("spec");

function code(f, p)
    return &+[Integers() | (Integers()!Coefficient(f, i))*p^i
        : i in [0..Degree(f)]];
end function;

function coordinates(P, p)
    if IsId(P) then return "O"; end if;
    return Join([Sprint(code(f, p)) : f in [Numerator(P[1]), Denominator(P[1]),
        Numerator(P[2]), Denominator(P[2])]], ":");
end function;

fprintf referencepath, "q,key,torsion,torsion_points,torsion_orders\n";
for p in [2, 3] do
    F<t> := FunctionField(GF(p));
    R<s, u> := PolynomialRing(GF(p), 2);
    weights := [2, 4, 6, 8, 12];
    equations := p eq 2 select [
        [F | 0, 0, 1, t^6, 0], [F | 0, 0, 1, t^7, 0],
        [F | t, 0, t^6+1, 0, 0], [F | t^2, 0, 1, t^8, t^6],
        [F | 1, 0, t^6, 1, t^6], [F | t^2, 0, 1, 0, 0]] else [
        [F | 0, 0, 0, t^8, t+1], [F | 0, t^3+t, 0, t^8+1, 0],
        [F | 0, t^4, 0, 0, 2], [F | 0, t^4+1, 0, t^4, 0],
        [F | 1, t^4, t, t^4+t, -t^2]];
    key := p eq 2 select 1346994688 else 75547326549;
    Append(~equations, [F | Evaluate(a, [t, F!1])
        : a in aInvariantsFromKey(R, key)]);
    for equation in equations do
        a := equation;
        // Keys omit a1 and a3 in characteristic 3, so normalize this chart first.
        if p eq 3 then
            a := [F | 0, a[2]+a[1]^2/4, 0, a[4]+a[1]*a[3]/2, a[5]+a[3]^2/4];
        end if;
        W := [R | &+[R | Coefficient(Numerator(a[i]), j)*s^j*u^(weights[i]-j)
            : j in [0..Degree(Numerator(a[i]))]] : i in [1..5]];
        for chart in [W, [R | Evaluate(f, [u, s]) : f in W]] do
            E := EllipticCurve([F | Evaluate(f, [t, F!1]) : f in chart]);
            factors, points := FibrationTorsion(E);
            T, mp := TorsionSubgroup(E);
            assert factors eq Invariants(T);
            assert Seqset(points) eq {mp(g) : g in T};
            labels := [coordinates(P, p) : P in points];
            orders := [Sprint(Order(P)) : P in points];
            fprintf referencepath, "%o,%o,%o,%o,%o\n", p, FibrationKey(chart),
                IsEmpty(factors) select "[]" else Join([Sprint(n) : n in factors], " "),
                Join(labels, ";"), Join(orders, " ");
        end for;
    end for;
end for;
quit;
