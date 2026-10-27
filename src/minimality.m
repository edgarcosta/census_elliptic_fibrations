import "forms.m": binary_basis, coefficients, form_from_vector, is_form_power;
import "models.m": valid_model, normalize_model, prime_field_key, ambient_space, ambient_vector;

function places(R, cusp)
    s := R.1;
    t := R.2;
    if Characteristic(R) eq 2 then
        result := [ R | t,s,s+t ];
        if not cusp then Append(~result,s^2+s*t+t^2); end if;
    else
        result := [ R | t,s,s+t,s-t ];
        if not cusp then
            result cat:= [ R | s^2+t^2,s^2+s*t-t^2,s^2-s*t-t^2 ];
        end if;
    end if;
    Sort(~result,func< x,y | TotalDegree(x) ne TotalDegree(y)
        select TotalDegree(x)-TotalDegree(y)
        else prime_field_key(coefficients(x,TotalDegree(x)),Characteristic(R))
            -prime_field_key(coefficients(y,TotalDegree(y)),Characteristic(R)) >);
    return result;
end function;

function residues(pi, j, d)
    // Degree-d forms modulo pi^j suffice for the weighted local conditions.
    R := Parent(pi);
    V := VectorSpace(BaseRing(R),d+1);
    U := sub< V | [ V!coefficients(pi^j*m,d)
        : m in binary_basis(R,d-j*TotalDegree(pi)) ] >;
    Q := sub< V | [ ReduceVector(U,v) : v in Basis(V) ] >;
    return [ R | form_from_vector(R,d,v) : v in Q ];
end function;

function local_conditions(forms, pi)
    R := Parent(pi);
    first,second := Explode(forms);
    conditions := [* *];
    if Characteristic(R) eq 2 then
        if not IsDivisibleBy(first,pi) then return conditions; end if;
        for q in residues(pi,1,2) do
            r := q*first+q^2;
            if not IsDivisibleBy(second+r*first,pi^3) then continue; end if;
            for v in residues(pi,3,6) do
                c4 := q*second+(v+r*q)*first+r^2;
                c6 := r^3+v*second+v^2+r*v*first;
                Append(~conditions,<r,c4,c6>);
            end for;
        end for;
    else
        if not IsDivisibleBy(first,pi^2) then return conditions; end if;
        for r in residues(pi,2,4) do
            if IsDivisibleBy(second-first*r,pi^4) then
                Append(~conditions,<r,R!0,second*r+first*r^2+r^3>);
            end if;
        end for;
    end if;
    return conditions;
end function;

function nonminimal_locus(R, forms)
    p := Characteristic(R);
    V := ambient_space(R);
    result := [* *];
    cusp := IsZero(forms[1]) and IsZero(forms[2]);
    for pi in places(R,cusp) do
        conditions := local_conditions(forms,pi);
        if #conditions eq 0 then continue; end if;
        e := TotalDegree(pi);
        U6 := [ R | pi^6*m : m in binary_basis(R,12-6*e) ];
        for condition in conditions do
            r,c4,c6 := Explode(condition);
            if p eq 2 then
                U4 := [ R | pi^4*m : m in binary_basis(R,8-4*e) ];
                U := sub< V | [ V!(ambient_vector([ R | u,r*u ])) : u in U4 ]
                    cat [ V!(ambient_vector([ R | 0,u ])) : u in U6 ] >;
                x := V!(ambient_vector([ R | c4,r*c4+c6 ]));
            else
                U := sub< V | [ V!ambient_vector([ R | u ]) : u in U6 ] >;
                x := V!ambient_vector([ R | -c6 ]);
            end if;
            Append(~result,<x,U>);
        end for;
    end for;
    return result;
end function;

function bad_place_certificate(W)
    W := normalize_model(W);
    R := Universe(W);
    p := Characteristic(R);
    forms := p eq 2 select [ W[1],W[3] ] else [ W[2],W[4] ];
    cusp := IsZero(forms[1]) and IsZero(forms[2]);
    if cusp and (p eq 2 select is_form_power(W[4],4) and is_form_power(W[5],2)
        else is_form_power(W[5],3)) then
        return false,<"all">;
    end if;
    bad := [ R | ];
    for pi in places(R,cusp) do
        for condition in local_conditions(forms,pi) do
            r,c4,c6 := Explode(condition);
            if p eq 2 then
                nonminimal := IsDivisibleBy(W[4]+c4,pi^4)
                    and IsDivisibleBy(W[5]+r*W[4]+c6,pi^6);
            else
                nonminimal := IsDivisibleBy(W[5]+c6,pi^6);
            end if;
            if nonminimal then Append(~bad,pi); break; end if;
        end for;
    end for;
    return IsEmpty(bad),<"finite",bad>;
end function;

intrinsic IsEverywhereMinimal(W::SeqEnum) -> BoolElt, Tup
{Whether the equation W is minimal at every closed place of P1, followed
by a bad-place certificate. The certificate is <"finite",places> for
all finitely many nonminimal places, or <"all"> for nonminimality at
every place. Nonregular models are included. See also FibrationProperties}
    require valid_model(W): "W must have five homogeneous coefficients of degrees 2,4,6,8,12 over one supported ring";
    return bad_place_certificate(W);
end intrinsic;
