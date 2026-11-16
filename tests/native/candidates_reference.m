SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/_native_candidate_samples.m": candidate_samples;

function wire(values)
    return IsEmpty(values) select "[]" else Join([Sprint(v) : v in values], " ");
end function;

file := Open(referencepath, "w");
delete file;
fprintf referencepath, "q,key,rank,artin_tate_product,sha_reg,admissible_sha,admissible_regulators,split_method\n";
for p in [2, 3] do
    F<t> := FunctionField(GF(p));
    R<s, u> := PolynomialRing(GF(p), 2);
    weights := [2, 4, 6, 8, 12];
    for sample in candidate_samples(F) do
        W := [R | &+[R | Coefficient(Numerator(sample[1][i]), j)*s^j*u^(weights[i]-j)
            : j in [0..Degree(Numerator(sample[1][i]))]] : i in [1..5]];
        for chart in [W, [R | Evaluate(f, [u, s]) : f in W]] do
            E := EllipticCurve([F | Evaluate(f, [t, F!1]) : f in chart]);
            data := FibrationInvariants(E);
            sha, regulators := FibrationCandidates(data, FibrationLocalData(E));
            assert sha eq sample[2];
            fprintf referencepath, "%o,%o,%o,%o,%o,%o,%o,%o\n", p, FibrationKey(chart),
                data`Rank, data`ArtinTateProduct, data`ShaRegulatorProduct,
                wire(sha), wire(regulators), #sha eq 1 select "constraints" else "unresolved";
        end for;
    end for;
end for;
quit;
