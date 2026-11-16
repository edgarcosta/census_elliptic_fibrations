SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/_native_torsion_samples.m": torsion_samples;
import "tests/_native_derived_fields.m": derived_fields;

procedure write_reference(path)
    file := Open(path, "w");
    delete file;
    fprintf path, "q,key,rank,geometric_rank,picard_number,geometric_picard_number,P2,zeta_sign,point_counts,brauer_height,artin_tate_product,ns_splitting_degree,j_degree,j_inseparability_exponent,is_isotrivial,sha_reg\n";
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        weights := [2, 4, 6, 8, 12];
        samples := torsion_samples(F);
        if p eq 3 then Append(~samples, <[F | 0, 0, 0, 1, t^7], []>); end if;
        for sample in samples do
            a := sample[1];
            if p eq 3 then
                a := [F | 0, a[2]+a[1]^2/4, 0, a[4]+a[1]*a[3]/2, a[5]+a[3]^2/4];
            end if;
            W := [R | &+[R | Coefficient(Numerator(a[i]), j)*s^j*u^(weights[i]-j)
                : j in [0..Degree(Numerator(a[i]))]] : i in [1..5]];
            for chart in [W, [R | Evaluate(f, [u, s]) : f in W]] do
                E := EllipticCurve([F | Evaluate(f, [t, F!1]) : f in chart]);
                result := FibrationInvariants(E);
                fields := [Sprint(p), Sprint(FibrationKey(chart))] cat derived_fields(result);
                fprintf path, "%o\n", Join(fields, ",");
            end for;
        end for;
    end for;
end procedure;

write_reference(referencepath);
quit;
