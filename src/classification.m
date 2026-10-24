import "forms.m": coefficients, form_from_vector, is_form_power;
import "models.m": valid_model, normalize_model, model_discriminant;
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
    V := VectorSpace(k,p eq 2 select 22 else 13);
    if first eq 0 and second eq 0 then
        S := [* <V!0,V> *];
        if p eq 2 then
            Z := sub< V | [ V.i : i in [1..9 by 2] ]
                cat [ V.i : i in [10..22 by 2] ] >;
        else
            Z := sub< V | [ V.i : i in [1..13 by 3] ] >;
        end if;
        return S,[* <V!0,Z> *];
    end if;
    if p eq 2 then
        rows := [];
        for z in Basis(V) do
            c := form_from_vector(R,8,[ z[i] : i in [1..9] ]);
            e := form_from_vector(R,12,[ z[i] : i in [10..22] ]);
            Append(~rows,coefficients(first^4*c^2+first^5*second*c+first^6*e,24));
        end for;
        linear := Matrix(k,22,25,rows);
        ok,x,U := IsConsistent(linear,Vector(k,coefficients(second^4+first^3*second^3,24)));
        S := ok select [* <V!x,U> *] else [* *];
    elif first ne 0 then
        ok,B := IsDivisibleBy(second^2*(first^2-second),first^3);
        S := ok select [* <V!coefficients(B,12),sub< V | >> *] else [* *];
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
