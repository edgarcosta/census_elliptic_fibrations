function derived_fields(result)
    return [Sprint(result`Rank), Sprint(result`GeometricRank),
        Sprint(result`PicardNumber), Sprint(result`GeometricPicardNumber),
        Join([Sprint(c) : c in Coefficients(result`P2)], " "),
        Sprint(result`ZetaSign), Join([Sprint(c) : c in result`PointCounts], " "),
        result`BrauerHeight eq Infinity() select "infinity" else Sprint(result`BrauerHeight),
        Sprint(result`ArtinTateProduct), Sprint(result`NeronSeveriSplittingDegree),
        Sprint(result`JDegree), Sprint(result`JInseparabilityExponent),
        result`IsIsotrivial select "true" else "false", Sprint(result`ShaRegulatorProduct)];
end function;
