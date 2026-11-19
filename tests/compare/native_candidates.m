import "tests/support/_support.m": compare_native;
import "tests/support/_native_candidate_samples.m": candidate_samples;

function wire(values)
    return IsEmpty(values) select "[]" else Join([Sprint(v) : v in values], " ");
end function;

referencepath, reference := TemporaryFile();
delete reference;
supplied := eval (Sprintf("referencepath := \"%o\";\n", referencepath)
    cat Read("tests/native_candidates_magma.m") cat "return Read(referencepath);");
Remove(referencepath);

procedure write_reference(path)
    fprintf path cat ".supplied", "%o", supplied;
    file := Open(path, "w");
    delete file;
    fprintf path, "q,key,admissible_sha,admissible_regulators\n";
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        for sample in candidate_samples(F) do
            W := aInvariants(EllipticCurve(sample[1]), R);
            for chart in [W, [R | Evaluate(f, [u, s]) : f in W]] do
                E := EllipticCurve(chart, F);
                data := FibrationInvariants(E);
                sha, regulators := FibrationCandidates(data, FibrationLocalData(E));
                assert sha eq sample[2];
                fprintf path, "%o,%o,%o,%o\n", p, FibrationKey(chart),
                    wire(sha), wire(regulators);
            end for;
        end for;
    end for;
end procedure;

compare_native(write_reference, "tests/compare/candidates.py" :
    Targets := "build/tests/native/candidates-exact", Suffixes := [".supplied"]);
