// Matches prefix orbits under base changes and quadratic shifts.
function reference_f2_prefix_match(R, old, new)
    pairs := [<Matrix(g),EF2FormFromMask(R,2,a)> : g in GL(2,GF(2)),a in [0..7]];
    canonical := function(P)
        return Min([EF2BitKey(EF2Coefficients(B[1],2) cat EF2Coefficients(B[3],6))
        where B is TransformNormalizedF2([R|P[1],0,P[2],0,0],h[1],h[2],R!0) : h in pairs]);
    end function;
    lookup := AssociativeArray();
    for i in [1..#old] do lookup[canonical(old[i])] := i; end for;
    match := [lookup[canonical(P)] : P in new];
    assert #new eq #old and #Keys(lookup) eq #old and #Seqset(match) eq #old;
    return match;
end function;
