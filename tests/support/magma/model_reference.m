function encode_fibres(rows)
    entries := [Strings() |];
    for r in rows do
        cycles := IsEmpty(r`ComponentCycles) select "[]"
            else Join([Sprint(n) : n in r`ComponentCycles], " ");
        fields := [Sprint(r`Place), Sprint(r`Degree), r`Type, Sprint(r`Index),
            Sprint(r`DiscriminantValuation), Sprint(r`Conductor),
            Sprint(r`Tamagawa), Sprint(r`Trace), cycles];
        Append(~entries, Join(fields, ":"));
    end for;
    return IsEmpty(entries) select "[]" else Join(entries, ";");
end function;

import "../_native_model_samples.m": model_samples;
import "../_native_fixed_models.m": fixed_models;

function coefficient_list(a, z)
    P := Parent(z);
    f := P!Evaluate(a, [z, P!1]);
    values := IsZero(f) select [BaseRing(P)!0] else Eltseq(f);
    return Join([IntegerToString(Integers()!c) : c in values], " ");
end function;

// The fixed forms pin all six cells; monomials pin every digit on both charts.
procedure write_model_reference(path)
    empty := Open(path, "w");
    delete empty;
    fprintf path, "q,key,kind,minimality,a1,a2,a3,a4,a6,infinity_a1,infinity_a2,infinity_a3,infinity_a4,infinity_a6,fibres,twist_fibres,degree,sign,height,twist_height,tamagawa_exponent\n";
    for p in [2, 3] do
        R<s, t> := PolynomialRing(GF(p), 2);
        P<z> := PolynomialRing(GF(p));
        degrees := [2, 4, 6, 8, 12];
        indices := p eq 2 select [1, 3, 4, 5] else [2, 4, 5];
        width := p eq 2 select 32 else 27;
        models := model_samples(R) cat fixed_models(R);
        models cat:= [[R | Evaluate(a, [t, s]) : a in W] : W in models];
        for key in [0, 1, p^width-1] do
            Append(~models, aInvariantsFromKey(R, key));
        end for;
        for i in indices do
            for j in [0..degrees[i]] do
                W := [R | 0, 0, 0, 0, 0];
                W[i] := (p-1)*s^(degrees[i]-j)*t^j;
                Append(~models, W);
            end for;
        end for;
        // Quadratic bad places require residue classes beyond rational points.
        for pi in (p eq 2 select [R | s^2+s*t+t^2] else
            [R | s^2+t^2, s^2+s*t-t^2, s^2-s*t-t^2]) do
            if p eq 2 then
                a1 := pi;
                a3 := pi^3;
                a4 := pi^4;
                a6 := pi^6;
                q := s*t;
                v := s^5*t+t^6;
                r := q*a1+q^2;
                W := [R | a1, 0, a3+r*a1,
                    a4+q*a3+(v+r*q)*a1+r^2,
                    a6+r*a4+r^3+v*a3+r*v*a1+v^2];
            else
                r := s^3*t+t^4;
                W := [R | 0, pi^2, 0, pi^4-r*pi^2,
                    pi^6+r*pi^4+r^2*pi^2+r^3];
            end if;
            assert not IsEverywhereMinimal(W);
            Append(~models, W);
        end for;
        // Full periods modulo 2^32 and 3^27 require multiplier 1 mod 12
        // and an increment coprime to 6.
        key := 17;
        for i in [1..100] do
            key := (1103515249*key+12347) mod p^width;
            Append(~models, aInvariantsFromKey(R, key));
        end for;
        for W in models do
            kind, minimality := FibrationProperties(W);
            row := [IntegerToString(p), IntegerToString(FibrationKey(W)), kind, minimality]
                cat [coefficient_list(a, z) : a in W]
                cat [coefficient_list(Evaluate(a, [t, s]), z) : a in W];
            if kind eq "smooth" and minimality eq "everywhere minimal" then
                F := FunctionField(GF(p));
                E := EllipticCurve(W, F);
                original := FibrationLocalData(E);
                twisted := FibrationLocalData(E : FrobeniusTwist := true);
                height := &+[r`Degree*r`DiscriminantValuation : r in original] div 12;
                twist_height := &+[r`Degree*r`DiscriminantValuation : r in twisted] div 12;
                exponent := &+[Valuation(r`Tamagawa, p) : r in twisted]
                    - &+[Valuation(r`Tamagawa, p) : r in original];
                sign := IsEven(twist_height-height+exponent) select 1 else -1;
                degree := &+[r`Degree*r`Conductor : r in original]-4;
                row cat:= [encode_fibres(original), encode_fibres(twisted), Sprint(degree),
                    Sprint(sign), Sprint(height), Sprint(twist_height), Sprint(exponent)];
            else
                row cat:= ["" : i in [1..7]];
            end if;
            fprintf path, "%o\n", Join(row, ",");
        end for;
    end for;
end procedure;
