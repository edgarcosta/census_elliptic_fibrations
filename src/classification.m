import "forms.m": coefficients, binary_basis, is_form_power;
import "models.m": valid_model, normalize_model, model_discriminant, ambient_space, ambient_vector, ambient_forms;
import "minimality.m": bad_place_certificate;

function generic_fibre_kind(W)
    W := normalize_model(W);
    if not IsZero(model_discriminant(W)) then return "smooth"; end if;
    if Characteristic(Universe(W)) eq 2 then
        cusp := W[1] eq 0 and W[3] eq 0;
        regular := not is_form_power(W[4],2) or not is_form_power(W[5],2);
    else
        cusp := W[2] eq 0 and W[4] eq 0;
        regular := not is_form_power(W[5],3);
    end if;
    return cusp and regular select "regular cusp" else "nonregular";
end function;

function generic_loci(R, forms)
    k := BaseRing(R);
    first,second := Explode(forms);
    p := Characteristic(R);
    V := ambient_space(R);
    if first eq 0 and second eq 0 then
        S := [* <V!0,V> *];
        if p eq 2 then
            Z := sub< V | [ V!ambient_vector([ R | m^2,0 ]) : m in binary_basis(R,4) ]
                cat [ V!ambient_vector([ R | 0,m^2 ]) : m in binary_basis(R,6) ] >;
        else
            Z := sub< V | [ V!ambient_vector([ R | m^3 ]) : m in binary_basis(R,4) ] >;
        end if;
        return S,[* <V!0,Z> *];
    end if;
    if p eq 2 then
        rows := [];
        for z in Basis(V) do
            c,e := Explode(ambient_forms(R,z));
            Append(~rows,coefficients(first^4*c^2+first^5*second*c+first^6*e,24));
        end for;
        linear := Matrix(k,Dimension(V),25,rows);
        ok,x,U := IsConsistent(linear,Vector(k,coefficients(second^4+first^3*second^3,24)));
        S := ok select [* <V!x,U> *] else [* *];
    elif first ne 0 then
        ok,B := IsDivisibleBy(second^2*(first^2-second),first^3);
        S := ok select [* <V!ambient_vector([ R | B ]),sub< V | >> *] else [* *];
    else
        S := [* *];
    end if;
    return S,S;
end function;

intrinsic FibrationProperties(W::SeqEnum) -> MonStgElt, MonStgElt
{The generic-fibre kind of W: "smooth", "regular cusp", or "nonregular",
followed by "everywhere minimal" or "nonminimal somewhere" for its equation.
Nonregular models are included. See also IsEverywhereMinimal, NumberOfFibrations}
    require valid_model(W): "W must have five homogeneous coefficients of degrees 2,4,6,8,12 over one supported ring";
    minimal,certificate := bad_place_certificate(W);
    return generic_fibre_kind(W),minimal select "everywhere minimal" else "nonminimal somewhere";
end intrinsic;
