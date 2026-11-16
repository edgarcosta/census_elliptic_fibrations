SetQuitOnError(true);
SetNthreads(1);
SetColumns(0);
AttachSpec("spec");
import "tests/_native_model_samples.m": model_samples;
import "tests/_native_fixed_models.m": fixed_models;
import "tests/_native_candidate_samples.m": candidate_samples;

file := Open(referencepath, "w");
delete file;
Write(referencepath, FibrationInvariantCSV());
for p in [2, 3] do
    R<s, u> := PolynomialRing(GF(p), 2);
    F<t> := FunctionField(GF(p));
    models := [sample[2] : sample in model_samples(R)] cat fixed_models(R);
    for sample in candidate_samples(F) do
        Append(~models, [R | &+[R | Coefficient(Numerator(sample[1][i]), j)
            *s^j*u^(2*[1,2,3,4,6][i]-j)
            : j in [0..Degree(Numerator(sample[1][i]))]] : i in [1..5]]);
    end for;
    for W in models do
        for chart in [W, [R | Evaluate(a, [u, s]) : a in W]] do
            Write(referencepath, FibrationInvariantCSV(FibrationInvariantRow(FibrationKey(chart), chart)));
        end for;
    end for;
end for;
quit;
