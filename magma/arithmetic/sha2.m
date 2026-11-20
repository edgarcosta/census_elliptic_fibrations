// Inputs come from a validated F_2 source row. Defect = -1 means that the
// alternating residue pairing has not supplied an exact defect.
function sha2_refinement(ordinary, descent_ok, f, v, rank, upper, orders : Defect := -1)
    if not ordinary then
        error if Defect ne -1, "Supersingular pairing is not applicable";
        return <"not_applicable", [], [], orders>;
    end if;
    if not descent_ok then
        error if Defect ne -1, "Pairing without successful saved descent";
        return <"not_computed", [], [], orders>;
    end if;
    error if not &and[Type(x) eq RngIntElt and x ge 0 : x in [f, v, rank]],
        "Invalid Selmer dimensions or Mordell-Weil rank";
    error if upper lt 0 or upper ne f + v - rank - 1, "Saved Selmer upper bound";

    defects := [];
    dimensions := {};
    survivors := {};
    for c in [0..Min(v, upper) by 2] do
        if Defect ne -1 and Defect ne c then
            continue;
        end if;
        d := f + v - rank - 1 - c;
        compatible := {n : n in orders | d eq 0 select IsOdd(n) else Valuation(n, 2) ge d};
        if not IsEmpty(compatible) then
            Append(~defects, c);
            Include(~dimensions, d);
            survivors join:= compatible;
        end if;
    end for;
    error if IsEmpty(defects), "Empty feasible defect intersection";
    return <#dimensions eq 1 select "exact" else "bound", Sort(Setseq(dimensions)),
        defects, [n : n in orders | n in survivors]>;
end function;
