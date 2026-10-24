for p in [2,3] do
    k := GF(p);
    R<s,t> := PolynomialRing(k,2);
    G := GL(2,k);
    for d in [2,4,6,8,12] do
        rho := EF_FormAction(R,G,d);
        V := VectorSpace(k,d+1);
        for g in G do
            for i in [1..d+1] do
                f := EF_FormFromVector(R,d,V.i*rho(g));
                assert f eq EF_BaseChange(s^(d-i+1)*t^(i-1),Matrix(g));
            end for;
        end for;
    end for;
end for;
