SetQuitOnError(true);
SetNthreads(1);

function forms(P, d)
    k := BaseRing(P);
    z := P.1;
    return [ P | &+[ P | k!v[i+1]*z^i : i in [0..d] ]
        : v in VectorSpace(k,d+1) ];
end function;

function coefficient_key(f, d)
    return &+[ Integers() | (Integers()!Coefficient(f,i))*3^i : i in [0..d] ];
end function;

function divisible(f, degree_bound, pi, exponent, infinity)
    if infinity then return Degree(f) le degree_bound-exponent; end if;
    return f mod pi^exponent eq 0;
end function;

function bad_equation(C,A,B,places,shifts)
    for i in [1..#places] do
        pi := places[i][1];
        infinity := places[i][2];
        if not divisible(C,4,pi,2,infinity) then continue; end if;
        for r in shifts[i] do
            if divisible(A-C*r,8,pi,4,infinity) and
                divisible(B+A*r+C*r^2+r^3,12,pi,6,infinity) then
                return true;
            end if;
        end for;
    end for;
    return false;
end function;

// Checks cusp, nonregular and bad-place coefficient sets independently.
procedure test_f3_coefficients()
    P<z> := PolynomialRing(GF(3));
    V4 := forms(P,4);
    V6 := forms(P,6);
    V8 := forms(P,8);
    places := [ <z,false>, <z-1,false>, <z+1,false>, <P!1,true>,
        <z^2+1,false>, <z^2+z-1,false>, <z^2-z-1,false> ];
    shifts := [ i eq 4 select [ r*z^3 : r in forms(P,1) ] else
        forms(P,2*Degree(places[i][1])-1) : i in [1..#places] ];
    cubes := { coefficient_key(r^3,12) : r in V4 };
    assert #cubes eq 3^5;

    cusp_bad_union := {};
    local_cusp_sets := [];
    for i in [1..4] do
        pi := places[i][1];
        local_set := { coefficient_key((i eq 4 select h else pi^6*h)-r^3,12)
            : h in V6, r in shifts[i] };
        assert #local_set eq 3^9;
        assert #(local_set meet cubes) eq 3^5;
        Include(~local_cusp_sets, local_set);
        cusp_bad_union join:= local_set;
    end for;
    for i in [1..3] do
        for j in [i+1..4] do
            assert local_cusp_sets[i] meet local_cusp_sets[j] eq cubes;
        end for;
    end for;
    assert #cusp_bad_union eq 4*3^9-3*3^5;
    cusp_count := 0;
    cusp_bad := 0;
    for b in [0..3^13-1] do
        if b in cubes then continue; end if;
        cusp_count +:= 1;
        if b in cusp_bad_union then cusp_bad +:= 1; end if;
    end for;
    assert cusp_count eq 3^13-3^5;
    assert cusp_bad eq 4*(3^9-3^5);

    nonregular := #cubes;
    nonregular_bad := #cubes;
    for C in V4 do
        if C eq 0 then continue; end if;
        C2 := C^2;
        C3 := C^3;
        for A in V8 do
            B, remainder := Quotrem(C2*A^2-A^3,C3);
            if remainder ne 0 or Degree(B) gt 12 then continue; end if;
            nonregular +:= 1;
            u, remainder := Quotrem(A,C);
            assert remainder eq 0 and Degree(u) le 4;
            assert B eq C*u^2-u^3;
            if bad_equation(C,A,B,places,shifts) then nonregular_bad +:= 1; end if;
        end for;
    end for;
    assert nonregular eq 3^10;
    assert nonregular_bad eq 99*3^5;

    divisors := [ places[i][1] : i in [5..7] ];
    for i in [1..3] do
        for j in [i+1..4] do
            Append(~divisors,places[i][1]*places[j][1]);
        end for;
    end for;
    zero_translates := { coefficient_key(r^3,12) : r in V4 };
    degree_two_sets := [];
    for D in divisors do
        pivot := Degree(D^2);
        exponents := [ i : i in [0..4] | i ne pivot ];
        residues := [ &+[ P | v[i]*z^exponents[i] : i in [1..4] ]
            : v in VectorSpace(GF(3),4) ];
        equations := {};
        singular_count := 0;
        for c,a,b in GF(3) do
            C := c*D^2;
            A0 := a*D^4;
            B0 := b*D^6;
            for r in residues do
                A := A0-C*r;
                B := B0+A0*r+C*r^2+r^3;
                Include(~equations,(coefficient_key(C,4)*3^9+
                    coefficient_key(A,8))*3^13+coefficient_key(B,12));
                if C^2*A^2-A^3-C^3*B eq 0 then
                    singular_count +:= 1;
                    assert C ne 0 or (A eq 0 and coefficient_key(B,12) in cubes);
                end if;
            end for;
        end for;
        assert #equations eq 3^7;
        assert singular_count eq 3^6;
        Append(~degree_two_sets,equations);
    end for;
    for i in [1..#degree_two_sets-1] do
        for j in [i+1..#degree_two_sets] do
            assert degree_two_sets[i] meet degree_two_sets[j] eq zero_translates;
        end for;
    end for;

    // Local translations are tested on jets, without global equation quotients.
    for C in [ z^2*f : f in forms(P,2) ] do
        signatures := { <coefficient_key((-C*r) mod z^4,3),
            coefficient_key((C*r^2+r^3) mod z^6,5)> : r in forms(P,1) };
        assert #signatures eq 9;
    end for;
end procedure;

try
    test_f3_coefficients();
catch e
    printf "MISMATCH: F3 coefficient sets %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
