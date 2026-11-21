import "../_native_torsion_samples.m": torsion_samples;
function derived_fields(result)
    return [Sprint(result`Rank), Sprint(result`GeometricRank),
        Sprint(result`PicardNumber), Sprint(result`GeometricPicardNumber),
        Join([Sprint(c) : c in Coefficients(result`P2)], " "),
        Sprint(result`ZetaSign), Join([Sprint(c) : c in result`PointCounts], " "),
        result`BrauerHeight eq Infinity() select "infinity" else Sprint(result`BrauerHeight),
        Sprint(result`ArtinTateProduct), Sprint(result`NeronSeveriSplittingDegree),
        Sprint(result`JDegree), Sprint(result`JInseparabilityExponent),
        Sprint(result`ShaRegulatorProduct)];
end function;

procedure write_reference(path)
    file := Open(path, "w");
    delete file;
    fprintf path, "q,key,rank,geometric_rank,picard_number,geometric_picard_number,P2,zeta_sign,point_counts,brauer_height,artin_tate_product,ns_splitting_degree,j_degree,j_inseparability_exponent,sha_reg\n";
    for p in [2, 3] do
        F<t> := FunctionField(GF(p));
        R<s, u> := PolynomialRing(GF(p), 2);
        samples := torsion_samples(F);
        if p eq 3 then Append(~samples, [F | 0, 0, 0, 1, t^7]); end if;
        for a in samples do
            W := aInvariantsFromKey(R, FibrationKey(aInvariants(EllipticCurve(a), R)));
            for chart in [W, [R | Evaluate(f, [u, s]) : f in W]] do
                E := EllipticCurve(chart, F);
                result := FibrationInvariants(E);
                fields := [Sprint(p), Sprint(FibrationKey(chart))] cat derived_fields(result);
                fprintf path, "%o\n", Join(fields, ",");
            end for;
        end for;
    end for;
end procedure;
