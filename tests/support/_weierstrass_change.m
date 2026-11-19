function homogeneous_forms(R, d)
    return [&+[v[j+1]*R.1^(d-j)*R.2^j : j in [0..d]]
        : v in VectorSpace(BaseRing(R), d+1)];
end function;

function equation_change(A, u, r, q, v)
    S<x, y> := PolynomialRing(Universe(A), 2);
    equation := y^2+A[1]*x*y+A[3]*y-x^3-A[2]*x^2-A[4]*x-A[5];
    changed := u^-6*Evaluate(equation, [u^2*x+r, u^3*y+q*u^2*x+v]);
    return [Universe(A) | MonomialCoefficient(changed, x*y),
        -MonomialCoefficient(changed, x^2), MonomialCoefficient(changed, y),
        -MonomialCoefficient(changed, x), -MonomialCoefficient(changed, S!1)];
end function;
