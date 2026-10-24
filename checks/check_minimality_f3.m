// Standalone validation, including an independent Magma Tate-algorithm check.
// Run from the repository root: magma -b checks/check_minimality_f3.m
load "k3_minimality_f3.m";

function K3F3_PolynomialValuation(f,pi)
    // Used only for nonzero discriminants and their nonzero numerators.
    assert f ne 0;
    n := 0;
    while f mod pi eq 0 do
        f := f div pi;
        n +:= 1;
    end while;
    return n;
end function;

function K3F3_HomogenizeTriple(R,W)
    assert #W eq 3;
    degrees := [4,8,12];
    s := R.1;
    t := R.2;
    assert &and[ Degree(W[i]) le degrees[i] : i in [1..3] ];
    return [R | &+[R | Coefficient(W[i],j)*s^j*t^(degrees[i]-j)
                       : j in [0..degrees[i]]] : i in [1..3]];
end function;

function K3F3_PrefixHasTranslationBrute(C,A,pp,residues)
    if C mod pp[2] ne 0 then return false; end if;
    return &or[ (A-C*r) mod pp[4] eq 0 : r in residues ];
end function;

function K3F3_CheckAgainstTate(W : Context := K3F3_default_context)
    // Checks all seven candidate places, including those with v(Delta)<12.
    // Return the number of checked places.  W must have Delta != 0.
    //
    // Documented API:
    //   info, Emin := LocalInformation(E, f);
    // f is an irreducible polynomial as a rational-function-field element,
    // or 1/T for infinity. info[2] is the local MINIMAL discriminant value.
    // https://magma.maths.usyd.edu.au/magma/handbook/text/1575
    //
    // Do not compare with Valuation(Delta(E),infinity) directly: E is
    // expressed in the affine chart, whereas the homogeneous O(2) model
    // has integral discriminant valuation 24-Degree(Delta(E)) there.
    R := K3F3_CheckModel(W);
    U := Context[1];
    u := U.1;
    F<T> := FunctionField(GF(3));
    affine := [U | Evaluate(a,[u,U!1]) : a in W];
    delta := K3F3_Discriminant(affine);
    assert delta ne 0;
    E := EllipticCurve([F | 0,affine[1],0,affine[2],affine[3]]);
    assert Discriminant(E) eq F!delta;
    ok,bad := IsEverywhereMinimalF3(W : Context := Context);
    prefix_auto := K3F3_PrefixAutomaticallyMinimal(W[1],W[2] :
                                                  Context := Context);
    if prefix_auto then assert ok; end if;
    any_prefix_candidate := false;
    tate_bad := [R | ];
    s := R.1;
    t := R.2;

    for i in [1..#Context[2]] do
        pi := Context[2][i];
        prefix_candidate := K3F3_PrefixHasTranslationBrute(
            affine[1],affine[2],Context[3][i],Context[4][Degree(pi)]);
        assert K3F3_PrefixCanBeNonMinimal(affine[1],affine[2],
                         Context[3][i]) eq prefix_candidate;
        any_prefix_candidate := any_prefix_candidate or prefix_candidate;
        info,Emin := LocalInformation(E,F!pi);
        oldv := K3F3_PolynomialValuation(delta,pi);
        newdelta := Discriminant(Emin);
        newv := K3F3_PolynomialValuation(U!Numerator(newdelta),pi)
              - K3F3_PolynomialValuation(U!Denominator(newdelta),pi);
        assert newv eq info[2];
        assert oldv ge newv and (oldv-newv) mod 12 eq 0;
        local_bad := K3F3_NonMinimalWithData(affine,Context[3][i],
                                            Context[4][Degree(pi)]);
        assert local_bad eq (oldv gt newv);
        if oldv gt newv then
            d := Degree(pi);
            Append(~tate_bad,&+[R | Coefficient(pi,j)*s^j*t^(d-j)
                                   : j in [0..d]]);
        end if;
    end for;

    info,Emin := LocalInformation(E,1/T);
    oldv := 24-Degree(delta);
    newdelta := Discriminant(Emin);
    newv := Degree(Denominator(newdelta))-Degree(Numerator(newdelta));
    assert newv eq info[2];
    assert oldv ge newv and (oldv-newv) mod 12 eq 0;
    infinity := [U | Evaluate(a,[U!1,u]) : a in W];
    prefix_candidate := K3F3_PrefixHasTranslationBrute(
        infinity[1],infinity[2],Context[3][1],Context[4][1]);
    assert K3F3_PrefixCanBeNonMinimal(infinity[1],infinity[2],
                     Context[3][1]) eq prefix_candidate;
    any_prefix_candidate := any_prefix_candidate or prefix_candidate;
    assert prefix_auto eq (not any_prefix_candidate);
    local_bad := K3F3_NonMinimalWithData(infinity,Context[3][1],
                                        Context[4][1]);
    assert local_bad eq (oldv gt newv);
    if oldv gt newv then Append(~tate_bad,t); end if;

    assert SequenceToSet(bad) eq SequenceToSet(tate_bad);
    assert ok eq IsEmpty(tate_bad);
    return #Context[2]+1;
end function;

procedure K3F3_RunMinimalityChecks(: RandomModels := 100, Seed := 31003)
    assert RandomModels ge 0;
    SetSeed(Seed);
    context := K3F3_Context();
    U := context[1];
    u := U.1;
    R<s,t> := PolynomialRing(GF(3),2);
    places_checked := 0;
    models_checked := 0;
    start := Cputime();

    // Good reduction model -> one deliberate pi scaling -> x translation.
    // Degree <= 2 guarantees that homogenization fits degrees 4,8,12.
    for pi in context[2] do
        d := Degree(pi);
        C := pi^2;
        for r in [U | 0,1,u^(2*d-1)+1] do
            local_model := [U | C,-C*r,pi^6+C*r^2+r^3];
            assert K3F3_NonMinimal(local_model,pi);
            W := K3F3_HomogenizeTriple(R,local_model);
            places_checked +:= K3F3_CheckAgainstTate(W : Context := context);
            models_checked +:= 1;
            // Swapping the homogeneous coordinates also exercises infinity.
            Wswap := [R | Evaluate(a,[t,s]) : a in W];
            places_checked +:= K3F3_CheckAgainstTate(Wswap : Context := context);
            models_checked +:= 1;
        end for;

        // Wild minimal example: Delta=-pi^12, but no x translation can
        // remove the pi term in B.  Invariants c4=c6=0 are insufficient.
        local_model := [U | 0,pi^4,pi];
        assert not K3F3_NonMinimal(local_model,pi);
        W := K3F3_HomogenizeTriple(R,local_model);
        places_checked +:= K3F3_CheckAgainstTate(W : Context := context);
        models_checked +:= 1;
    end for;

    // Singular generic fibre is rejected before EllipticCurve is called.
    ok,bad := IsEverywhereMinimalF3([R | 0,0,s^12] : Context := context);
    assert not ok and IsEmpty(bad);

    degrees := [4,8,12];
    for i in [1..RandomModels] do
        repeat
            W := [R | &+[R | Random(GF(3))*s^j*t^(d-j) : j in [0..d]]
                      : d in degrees];
        until K3F3_Discriminant(W) ne 0;
        places_checked +:= K3F3_CheckAgainstTate(W : Context := context);
        models_checked +:= 1;
    end for;
    printf "Minimality/Tate checks passed: %o models, %o local checks, %o CPU s\n",
           models_checked,places_checked,Cputime(start);
end procedure;

K3F3_RunMinimalityChecks();
quit;
