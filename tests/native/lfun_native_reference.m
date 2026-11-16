SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/_native_fixed_models.m": fixed_models;

procedure write_lfun_reference(path)
    fprintf path, "q,key,degree,sign,lower_L,L_polynomial,trace_sums\n";
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        models := fixed_models(R);
        key := 17;
        for i in [1..8] do
            key := (1103515249*key+12347) mod p^(p eq 2 select 32 else 27);
            W := aInvariantsFromKey(R, key);
            kind, minimality := FibrationProperties(W);
            if kind eq "smooth" and minimality eq "everywhere minimal" then
                Append(~models, W);
            end if;
        end for;
        for W in models do
            for chart in [W, [R | Evaluate(a, [u, s]) : a in W]] do
                E := EllipticCurve([F | Evaluate(a, [t, F!1]) : a in chart]);
                result := FibrationLPolynomial(E);
                lists := [result`Lower, Coefficients(result`Polynomial), result`TraceSums];
                encoded := [IsEmpty(values) select "[]"
                    else Join([Sprint(c) : c in values], " ") : values in lists];
                fprintf path, "%o,%o,%o,%o,%o\n", p, FibrationKey(chart),
                    result`Degree, result`Sign, Join(encoded, ",");
            end for;
        end for;
    end for;
end procedure;

write_lfun_reference(referencepath);
quit;
