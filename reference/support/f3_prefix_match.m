// Matches prefix orbits under base changes and quartic shifts in A.
function reference_f3_prefix_match(R, old, new)
    k := BaseRing(R);
    V8 := VectorSpace(k,9);
    B4 := EF_BinaryBasis(R,4);
    base := EF3_BaseMatrices(k);
    function canonical_prefix(C,A)
        candidates := [];
        for M in base do
            D := EF_BaseChange(C,M);
            S := sub< V8 | [ V8!EF_Coefficients(D*r,8) : r in B4 ] >;
            a := ReduceVector(S,V8!EF_Coefficients(EF_BaseChange(A,M),8));
            Append(~candidates,<EF3_Key(EF_Coefficients(D,4)),EF3_Key(a)>);
        end for;
        return Min(candidates);
    end function;
    index := AssociativeArray();
    for i in [1..#old] do
        index[canonical_prefix(old[i][1],old[i][2])] := i;
    end for;
    assert #Keys(index) eq #old;
    match := [ index[canonical_prefix(P[1],P[2])] : P in new ];
    assert #new eq #old and #SequenceToSet(match) eq #old;
    return match;
end function;
