import "../_native_model_samples.m": model_samples;
import "../_native_fixed_models.m": fixed_models;
import "../_native_candidate_samples.m": candidate_samples;

procedure write_records_reference(path)
    file := Open(path, "w");
    delete file;
    Write(path, FibrationInvariantCSV());
    for p in [2, 3] do
        R<s, u> := PolynomialRing(GF(p), 2);
        F<t> := FunctionField(GF(p));
        models := model_samples(R) cat fixed_models(R);
        for sample in candidate_samples(F) do
            Append(~models, aInvariants(EllipticCurve(sample[1]), R));
        end for;
        for W in models do
            for chart in [W, [R | Evaluate(a, [u, s]) : a in W]] do
                Write(path, FibrationInvariantCSV(FibrationInvariantRow(FibrationKey(chart), chart)));
            end for;
        end for;
    end for;
end procedure;
