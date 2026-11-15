import "../src/models.m": prime_field_key;

LOCAL_DATA_FORMAT := recformat<Place, Degree, Type, Index,
    DiscriminantValuation, Conductor, Tamagawa, Trace, ComponentCycles>;

function component_cycles(li)
    M := FrobeniusActionOnReducibleFiber(li);
    used := {};
    cycles := [Integers() |];
    for i in [1..Nrows(M)] do
        if i in used then continue; end if;
        j := i;
        length := 0;
        repeat
            Include(~used, j);
            length +:= 1;
            j := Index(Eltseq(M[j]), 1);
        until j eq i;
        Append(~cycles, length);
    end for;
    Sort(~cycles);
    return cycles;
end function;

intrinsic FibrationLocalData(E::CrvEll[FldFunRat] : FrobeniusTwist := false) -> SeqEnum
{The bad-place fibre records of E over GF(2)(t) or GF(3)(t).
The parameter FrobeniusTwist (default false) uses coefficients a_i^p before
local minimization. See also LocalInformation}
    F := BaseRing(E);
    k := BaseRing(F);
    require Type(k) eq FldFin and #k in [2, 3]:
        "constant field must be GF(2) or GF(3)";
    if FrobeniusTwist then
        E := EllipticCurve([a^Characteristic(k) : a in aInvariants(E)]);
    end if;
    P := PolynomialRing(k);
    rows := [];
    for li in LocalInformation(E) do
        infinity := not IsFinite(li[1]);
        if infinity then
            place := "inf";
            degree := 1;
        else
            f := P!UniformizingElement(li[1]);
            f /:= LeadingCoefficient(f);
            place := prime_field_key(Reverse(Coefficients(f)), #k);
            degree := Degree(f);
        end if;
        symbol := Sprint(li[5]);
        family := symbol;
        n := 0;
        if #symbol ge 2 and symbol[1] eq "I" and symbol[2] in "0123456789" then
            starred := symbol[#symbol] eq "*";
            family := starred select "I*" else "I";
            n := StringToInteger(symbol[2..#symbol-(starred select 1 else 0)]);
        end if;
        trace := family eq "I" select (li[6] select 1 else -1) else 0;
        Append(~rows, rec<LOCAL_DATA_FORMAT | Place := place, Degree := degree,
            Type := family, Index := n, DiscriminantValuation := li[2],
            Conductor := li[3], Tamagawa := li[4], Trace := trace,
            ComponentCycles := component_cycles(li)>);
    end for;
    Sort(~rows, func<a, b | a`Place cmpeq b`Place select 0
        else a`Place cmpeq "inf" select 1
        else b`Place cmpeq "inf" select -1 else a`Place-b`Place>);
    return rows;
end intrinsic;
