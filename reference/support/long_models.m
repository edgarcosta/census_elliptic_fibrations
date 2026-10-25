// Independent long-model formulas keep these checks separate from scan orbits.
function math_monomials(R,d) return [R| R.1^(d-i)*R.2^i : i in [0..d]]; end function;
function math_coefficients(f,d) return [MonomialCoefficient(f,m) : m in math_monomials(Parent(f),d)]; end function;
function math_form(R,v,d) B := math_monomials(R,d); return &+[R| v[i]*B[i] : i in [1..d+1]]; end function;
function math_base_change(f,M) R := Parent(f); return Evaluate(f,[M[1,1]*R.1+M[1,2]*R.2, M[2,1]*R.1+M[2,2]*R.2]); end function;
function math_forms(R,d) return [ math_form(R,v,d) : v in VectorSpace(BaseRing(R),d+1) ]; end function;
function math_base_matrices(k)
    mats := [ Matrix(k,2,2,[a,b,c,d]) : a,b,c,d in k ];
    mats := [ M : M in mats | Determinant(M) ne 0 ];
    return [ M : M in mats | [x : x in Eltseq(M) | x ne 0][1] eq 1 ];
end function;
function math_weierstrass_change(A,u,r,s,t)
    a1,a2,a3,a4,a6 := Explode(A);
    return [ Universe(A) | (a1+2*s)/u, (a2-s*a1+3*r-s^2)/u^2, (a3+r*a1+2*t)/u^3,
        (a4-s*a3+2*r*a2-(t+r*s)*a1+3*r^2-2*s*t)/u^4,
        (a6+r*a4+r^2*a2+r^3-t*a3-r*t*a1-t^2)/u^6 ];
end function;
function math_discriminant(A)
    a1,a2,a3,a4,a6 := Explode(A);
    b2 := a1^2+4*a2; b4 := a1*a3+2*a4; b6 := a3^2+4*a6;
    b8 := a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
    return -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6;
end function;
function math_places(R,maxdeg)
    k := BaseRing(R); res := [R|];
    for d in [1..maxdeg] do
        for f in math_forms(R,d) do
            if f eq 0 then continue; end if;
            c := math_coefficients(f,d); lead := [x : x in c | x ne 0][1];
            if lead ne 1 then continue; end if;
            if IsIrreducible(f) then Append(~res,f); end if;
        end for;
    end for;
    key := func< f | &+[ Integers()!c[i]*(#k)^(#c-i) : i in [1..#c] ] where c := math_coefficients(f,TotalDegree(f)) >;
    Sort(~res, func< f,g | TotalDegree(f) ne TotalDegree(g) select TotalDegree(f)-TotalDegree(g) else key(f)-key(g) >);
    return res;
end function;
// For places of degree at most two, these spaces cover q mod pi,
// r mod pi^2 and t mod pi^3. Higher lifts preserve weighted divisibility.
function math_local_nonminimal(A,pi,V4,V2,V6)
    a1,a2,a3,a4,a6 := Explode(A);
    if Characteristic(Universe(A)) eq 2 then
        if not IsDivisibleBy(a1,pi) then return false; end if;
        for q in V2 do
            r := a2+q*a1+q^2;
            if not IsDivisibleBy(a3+r*a1,pi^3) then continue; end if;
            for t in V6 do
                B := math_weierstrass_change(A,1,r,q,t);
                assert B[1] eq a1 and B[2] eq 0 and B[3] eq a3+r*a1;
                if IsDivisibleBy(B[4],pi^4) and IsDivisibleBy(B[5],pi^6) then return true; end if;
            end for;
        end for;
        return false;
    end if;
    s := a1;
    for r in V4 do
        t := a3+r*a1;
        B := math_weierstrass_change(A,1,r,s,t);
        assert B[1] eq 0 and B[3] eq 0;
        if IsDivisibleBy(B[2],pi^2) and IsDivisibleBy(B[4],pi^4) and IsDivisibleBy(B[5],pi^6) then return true; end if;
    end for;
    return false;
end function;
function math_bad_places(A,places,V4,V2,V6)
    return [ pi : pi in places | math_local_nonminimal(A,pi,V4,V2,V6) ];
end function;
