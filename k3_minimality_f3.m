// Local minimality for completed-square models over GF(3):
// y^2 = x^3 + C*x^2 + A*x + B, with homogeneous degrees 4, 8, 12.
// This file is self-contained.  It does not call Tate's algorithm.
// checks/check_minimality_f3.m compares it with Tate's algorithm.

function K3F3_CheckModel(W)
    assert #W eq 3;
    R := Parent(W[1]);
    assert Type(R) eq RngMPol and Rank(R) eq 2;
    assert Type(BaseRing(R)) eq FldFin and #BaseRing(R) eq 3;
    degrees := [4,8,12];
    for i in [1..3] do
        assert Parent(W[i]) eq R;
        for m in Monomials(W[i]) do
            assert TotalDegree(m) eq degrees[i];
        end for;
    end for;
    return R;
end function;

function K3F3_Discriminant(W)
    // Valid for homogeneous or univariate completed-square triples.
    C,A,B := Explode(W);
    return C^2*A^2 - A^3 - C^3*B;
end function;

function K3F3_Context()
    // <univariate ring, six finite places, powers, r representatives>.
    // Entry 4[d] is a list of representatives modulo a degree-d pi^2;
    // every polynomial of degree < 2*d occurs once.
    U<u> := PolynomialRing(GF(3));
    primes := [U | u,u-1,u-2,u^2+1,u^2+u+2,u^2+2*u+2];
    assert &and[ IsIrreducible(pi) : pi in primes ];
    powers := [[pi^j : j in [1..12]] : pi in primes];
    residues := [ [U |
        &+[U | (GF(3)!b[j+1])*u^j : j in [0..2*d-1]]
        : b in CartesianPower([0,1,2],2*d)] : d in [1,2] ];
    return <U,primes,powers,residues>;
end function;

K3F3_default_context := K3F3_Context();

function K3F3_NonMinimalWithData(W, pp, residues)
    C,A,B := Explode(W);
    if C mod pp[2] ne 0 then return false; end if;
    for r in residues do
        if (A-C*r) mod pp[4] ne 0 then continue; end if;
        if (B+A*r+C*r^2+r^3) mod pp[6] eq 0 then
            return true;
        end if;
    end for;
    return false;
end function;

function K3F3_PrefixCanBeNonMinimal(C, A, pp)
    // Does some r modulo pi^2 make C and A-C*r divisible by pi^2,pi^4?
    // Once v(C)>=2, this is equivalent to v(A)>=Min(v(C),4).
    // Testing divisibilities avoids computing a valuation of C=0.
    if C mod pp[2] ne 0 then return false; end if;
    if C mod pp[3] ne 0 then
        modulus := pp[2];
    elif C mod pp[4] ne 0 then
        modulus := pp[3];
    else
        modulus := pp[4];
    end if;
    return A mod modulus eq 0;
end function;

function K3F3_PrefixAutomaticallyMinimal(C, A :
                              Context := K3F3_default_context)
    // C,A are homogeneous of degrees 4,8 over GF(3)[s,t].
    // True iff no degree <= 2 place admits r modulo pi^2 satisfying
    // pi^2|C and pi^4|(A-C*r).  Then every generically smooth B-fibre
    // is everywhere minimal.  A false result only requests the B test.
    // This function does not test Delta != 0.
    R := Parent(C);
    R := K3F3_CheckModel([R | C,A,0]);
    U := Context[1];
    u := U.1;
    c := U!Evaluate(C,[u,U!1]);
    a := U!Evaluate(A,[u,U!1]);
    for i in [1..#Context[2]] do
        if K3F3_PrefixCanBeNonMinimal(c,a,Context[3][i]) then
            return false;
        end if;
    end for;
    c := U!Evaluate(C,[U!1,u]);
    a := U!Evaluate(A,[U!1,u]);
    return not K3F3_PrefixCanBeNonMinimal(c,a,Context[3][1]);
end function;

function K3F3_NonMinimal(W, pi)
    // W is an integral univariate triple and pi is monic irreducible,
    // of degree 1 or 2.  For Delta != 0, returns true exactly when the
    // given integral equation is nonminimal at pi.
    //
    // With x=X+r, the new coefficients are
    // C, A-C*r, B+A*r+C*r^2+r^3.  They must be divisible by
    // pi^2, pi^4, pi^6.  Once the first two conditions hold, replacing
    // r by r+pi^2*h does not affect the third condition modulo pi^6.
    assert #W eq 3;
    U := Parent(pi);
    assert Type(U) eq RngUPol;
    assert Type(BaseRing(U)) eq FldFin and #BaseRing(U) eq 3;
    assert &and[ Parent(a) eq U : a in W ];
    d := Degree(pi);
    assert d in [1,2];
    assert LeadingCoefficient(pi) eq 1 and IsIrreducible(pi);
    u := U.1;
    pp := [pi^j : j in [1..12]];
    residues := [U |
        &+[U | (BaseRing(U)!b[j+1])*u^j : j in [0..2*d-1]]
        : b in CartesianPower([0,1,2],2*d)];
    return K3F3_NonMinimalWithData(W,pp,residues);
end function;

function IsEverywhereMinimalF3(W : Context := K3F3_default_context)
    // Input: [C,A,B] over GF(3)[s,t], homogeneous of degrees 4,8,12.
    // Returns two values: Boolean, sequence of homogeneous bad-place forms.
    // Delta = 0 returns false and an empty sequence, as in the F2 code.
    // The Boolean is true precisely for an everywhere-minimal model
    // with smooth generic fibre (hence an elliptic K3 surface).
    R := K3F3_CheckModel(W);
    s := R.1;
    t := R.2;
    U := Context[1];
    u := U.1;
    affine := [U | Evaluate(a,[u,U!1]) : a in W];
    delta := K3F3_Discriminant(affine);
    bad := [R | ];
    if delta eq 0 then return false,bad; end if;

    // For a nonzero homogeneous Delta of degree 24, nonminimality
    // implies v_pi(Delta) >= 12, so only places of degree <= 2 occur.
    for i in [1..#Context[2]] do
        pi := Context[2][i];
        pp := Context[3][i];
        if delta mod pp[12] ne 0 then continue; end if;
        d := Degree(pi);
        if K3F3_NonMinimalWithData(affine,pp,Context[4][d]) then
            Append(~bad,&+[R | Coefficient(pi,j)*s^j*t^(d-j)
                              : j in [0..d]]);
        end if;
    end for;

    // At infinity, use the O(2) trivialization with (s,t)=(1,u).
    // Its discriminant valuation is 24-Degree(delta).
    if Degree(delta) le 12 then
        infinity := [U | Evaluate(a,[U!1,u]) : a in W];
        if K3F3_NonMinimalWithData(infinity,Context[3][1],
                                 Context[4][1]) then
            Append(~bad,t);
        end if;
    end if;
    return IsEmpty(bad),bad;
end function;
