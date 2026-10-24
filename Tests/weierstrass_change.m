for p in [2,3] do
    k := GF(p);
    F<T> := FunctionField(k);
    A := [ F | T+1,T^2+T,T^3+1,T^4+T,T^6+1 ];
    E := EllipticCurve(A);
    for u in [ F | 1,T ] do
        for r,s,t in [ F | 0,1,T ] do
            B := EF_WeierstrassChange(A,u,r,s,t);
            // Transformation maps old coordinates to new ones; the routine
            // substitutes old coordinates in terms of the new ones.
            Enew := Transformation(E,[-r/u^2,-s/u,(r*s-t)/u^3,1/u]);
            assert B eq aInvariants(Enew);
        end for;
    end for;

    R<z,w> := PolynomialRing(k,2);
    A := [ R | z^2,0,z^6,z^8,w^12 ];
    B := EF_WeierstrassChange(A,k!1,z^4,z^2,w^6);
    assert Universe(B) eq R;
    assert EF_CheckForms(B,[2,4,6,8,12]) eq R;

    P<x,y> := PolynomialRing(F,2);
    for A in [[ F | 0,0,0,0,0 ],aInvariants(E)] do
        u := T; r := T+1; s := T^2; t := T^3+1;
        B := EF_WeierstrassChange(A,u,r,s,t);
        oldx := u^2*x+r;
        oldy := u^3*y+s*u^2*x+t;
        oldeq := oldy^2+A[1]*oldx*oldy+A[3]*oldy
                 -oldx^3-A[2]*oldx^2-A[4]*oldx-A[5];
        neweq := y^2+B[1]*x*y+B[3]*y-x^3-B[2]*x^2-B[4]*x-B[5];
        assert oldeq eq u^6*neweq;
    end for;
end for;
