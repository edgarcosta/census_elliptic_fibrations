SetQuitOnError(true);
AttachSpec("spec");

function coefficient_forms(P,d)
    return [ P!Intseq(n,2) : n in [0..2^(d+1)-1] ];
end function;

function coefficient_discriminant(a1,a3,a4,a6)
    return a3^4+a1^3*a3^3+a1^4*a4^2+a1^5*a3*a4+a1^6*a6;
end function;

// Enumerates the complete nonminimal union and solves the singular equations.
procedure test_f2_coefficients()
    P<z> := PolynomialRing(GF(2));
    spaces := [ coefficient_forms(P,d) : d in [2,6,8,12] ];
    keys := [];
    for space in spaces do
        lookup := AssociativeArray();
        for i in [1..#space] do lookup[space[i]] := i-1; end for;
        Append(~keys,lookup);
    end for;
    function equation_key(a1,a3,a4,a6)
        return ((keys[1][a1]*128+keys[2][a3])*512+keys[3][a4])*8192+keys[4][a6];
    end function;
    function equation_from_key(key)
        a6 := spaces[4][1+key mod 8192]; key div:= 8192;
        a4 := spaces[3][1+key mod 512]; key div:= 512;
        a3 := spaces[2][1+key mod 128]; key div:= 128;
        return spaces[1][1+key],a3,a4,a6;
    end function;
    bad := { Integers() | };
    // Infinity lowers each coefficient's degree bound in the other affine chart.
    for place in [1..4] do
        infinity := place eq 3;
        d := place eq 4 select 2 else 1;
        pi := [ z,z+1,P!1,z^2+z+1 ][place];
        local_set := { Integers() | };
        for h1 in coefficient_forms(P,2-d), h3 in coefficient_forms(P,6-3*d),
            h4 in coefficient_forms(P,8-4*d), h6 in coefficient_forms(P,12-6*d) do
            b1 := pi*h1; b3 := pi^3*h3; b4 := pi^4*h4; b6 := pi^6*h6;
            for q0 in coefficient_forms(P,d-1), v0 in coefficient_forms(P,3*d-1) do
                q := infinity select q0*z^2 else q0;
                v := infinity select v0*z^4 else v0;
                r := q*b1+q^2;
                a3 := b3+r*b1;
                a4 := b4+q*b3+(v+r*q)*b1+r^2;
                a6 := b6+r*b4+r^3+v*b3+r*v*b1+v^2;
                Include(~local_set,equation_key(b1,a3,a4,a6));
            end for;
        end for;
        assert #local_set eq (d eq 1 select 2^22 else 2^12);
        bad join:= local_set;
    end for;
    assert #bad eq 3*2^22-2*2^12;
    bad_cells := [0,0,0];
    for key in bad do
        a1,a3,a4,a6 := equation_from_key(key);
        delta := coefficient_discriminant(a1,a3,a4,a6);
        kind := delta ne 0 select 1 else
            (a1 eq 0 and a3 eq 0 and (Derivative(a4) ne 0 or Derivative(a6) ne 0)
                select 2 else 3);
        bad_cells[kind] +:= 1;
    end for;
    assert bad_cells eq [12363776,190464,20480];
    nonregular := { Integers() | };
    for a1 in spaces[1] do
        if a1 eq 0 then continue; end if;
        for a3 in spaces[2], a4 in spaces[3] do
            a6,remainder := Quotrem(coefficient_discriminant(a1,a3,a4,P!0),a1^6);
            if remainder eq 0 and Degree(a6) le 12 then
                Include(~nonregular,equation_key(a1,a3,a4,a6));
            end if;
        end for;
    end for;
    for f in coefficient_forms(P,4), g in coefficient_forms(P,6) do
        Include(~nonregular,equation_key(P!0,P!0,f^2,g^2));
    end for;
    assert #nonregular eq 2^15 and #(nonregular meet bad) eq 20480;
    cusp := 0; cusp_bad := 0;
    for a4 in spaces[3], a6 in spaces[4] do
        if Derivative(a4) eq 0 and Derivative(a6) eq 0 then continue; end if;
        cusp +:= 1;
        if equation_key(P!0,P!0,a4,a6) in bad then cusp_bad +:= 1; end if;
    end for;
    assert cusp eq 2^22-2^12 and cusp_bad eq 190464;
    R := PolynomialRing(GF(2),2);
    n,cells,rows := NumberOfEquations(R);
    assert n eq 2^32;
    for i in [1..3] do
        kind := ["smooth","regular cusp","nonregular"][i];
        assert cells[<kind,"nonminimal somewhere">] eq bad_cells[i];
    end for;
    assert cells[<"regular cusp","everywhere minimal">] eq cusp-cusp_bad;
    assert cells[<"nonregular","everywhere minimal">] eq #nonregular-20480;
end procedure;

try
    test_f2_coefficients();
catch e
    printf "MISMATCH: F2 coefficient sets %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
