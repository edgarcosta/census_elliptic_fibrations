import "models.m": valid_ring;
import "f2.m": f2_prefixes, f2_fibre;
import "f3.m": f3_prefixes, f3_fibre;
import "engine.m": fibre_counts, fibre_representatives;

KINDS := [ "smooth","regular cusp","nonregular" ];
MINIMALITIES := [ "everywhere minimal","nonminimal somewhere" ];

function selection_error(indices, kinds, minimalities, n)
    if Type(indices) ne SeqEnum or Type(kinds) ne SeqEnum or Type(minimalities) ne SeqEnum then
        return "Selections must be sequences";
    end if;
    if not &and[ Type(i) eq RngIntElt and i ge 1 and i le n : i in indices ] then
        return "PrefixIndices must contain valid prefix indices";
    end if;
    if #Seqset(indices) ne #indices then return "PrefixIndices must be distinct"; end if;
    if not &and[ Type(k) eq MonStgElt and k in KINDS : k in kinds ] or #Seqset(kinds) ne #kinds then
        return "GenericFibreKinds must contain distinct values: smooth, regular cusp, nonregular";
    end if;
    if not &and[ Type(m) eq MonStgElt and m in MINIMALITIES : m in minimalities ] or #Seqset(minimalities) ne #minimalities then
        return "Minimalities must contain distinct values: everywhere minimal, nonminimal somewhere";
    end if;
    return "";
end function;

function selections(indices, kinds, minimalities, n)
    indices := #indices eq 0 select [1..n] else indices;
    kinds := #kinds eq 0 select KINDS else [ k : k in KINDS | k in kinds ];
    minimalities := #minimalities eq 0 select MINIMALITIES
        else [ m : m in MINIMALITIES | m in minimalities ];
    return indices,kinds,minimalities;
end function;

function prefixes(R)
    return Characteristic(R) eq 2 select f2_prefixes(R) else f3_prefixes(R);
end function;

function fibre(R, P)
    return Characteristic(R) eq 2 select f2_fibre(P) else f3_fibre(P);
end function;

function pattern_to_cell(pattern)
    assert pattern[1] or not pattern[2];
    kind := not pattern[1] select KINDS[1] else (pattern[2] select KINDS[3] else KINDS[2]);
    return <kind,pattern[3] select MINIMALITIES[2] else MINIMALITIES[1]>;
end function;

function empty_cells(kinds, minimalities)
    cells := AssociativeArray();
    for k in kinds do
        for m in minimalities do cells[<k,m>] := 0; end for;
    end for;
    return cells;
end function;

function selected_patterns(kinds, minimalities)
    return [ <s,r,n> : s,r,n in [false,true] | (s or not r)
        and pattern_to_cell(<s,r,n>)[1] in kinds
        and pattern_to_cell(<s,r,n>)[2] in minimalities ];
end function;

function selected_representative_keys(R, P, kinds, minimalities)
    F := fibre(R,P);
    tagged := fibre_representatives(F,selected_patterns(kinds,minimalities));
    cells := empty_cells(kinds,minimalities);
    keys := [ Integers() | ];
    for pair in tagged do
        Append(~keys,pair[1]);
        cells[pattern_to_cell(pair[2])] +:= 1;
    end for;
    return keys,cells;
end function;

intrinsic FibrationPrefixes(R::RngMPol) -> SeqEnum
{The ordered prefix pairs for the fibration census over R, a binary
polynomial ring over GF(2) or GF(3). The pairs are <a1,a3> over GF(2)
and <C,A> over GF(3), after completing the square. See also NumberOfFibrations}
    require valid_ring(R): "R must be a binary polynomial ring over GF(2) or GF(3)";
    return [ <P[1],P[2]> : P in prefixes(R) ];
end intrinsic;

intrinsic NumberOfFibrations(R::RngMPol : PrefixIndices := [],
    GenericFibreKinds := [], Minimalities := []) -> RngIntElt, Assoc, SeqEnum
{The number of selected model orbits over R, followed by their generic-fibre
and minimality cell counts and pairs <prefix index, cell counts>.
The parameter PrefixIndices (default []) selects all prefixes when empty.
The parameter GenericFibreKinds (default []) selects "smooth", "regular cusp",
and "nonregular" when empty, or a subset of these case-sensitive values.
The parameter Minimalities (default []) selects "everywhere minimal" and
"nonminimal somewhere" when empty, or a subset of these case-sensitive values.
See also FibrationProperties, FibrationRepresentativeKeys}
    require valid_ring(R): "R must be a binary polynomial ring over GF(2) or GF(3)";
    P := prefixes(R);
    message := selection_error(PrefixIndices,GenericFibreKinds,Minimalities,#P);
    require message eq "": message;
    indices,kinds,minimalities := selections(PrefixIndices,GenericFibreKinds,Minimalities,#P);
    totals := empty_cells(kinds,minimalities);
    rows := [];
    for i in indices do
        timestamp := Cputime();
        F := fibre(R,P[i]);
        patterns := fibre_counts(F);
        cells := empty_cells(kinds,minimalities);
        for pattern in Keys(patterns) do
            if not pattern[1] and pattern[2] then assert patterns[pattern] eq 0; continue; end if;
            cell := pattern_to_cell(pattern);
            if IsDefined(cells,cell) then cells[cell] +:= patterns[pattern]; end if;
        end for;
        for cell in Keys(totals) do totals[cell] +:= cells[cell]; end for;
        Append(~rows,<i,cells>);
        vprintf FibrationCensus,1: "Completed prefix %o\n",i;
        vprintf FibrationCensus,2: "Dimension %o, group order %o, CPU %o s\n",
            F`Dimension,#F`Group,Cputime(timestamp);
    end for;
    return &+[ Integers() | totals[cell] : cell in Keys(totals) ],totals,rows;
end intrinsic;

intrinsic FibrationRepresentativeKeys(R::RngMPol, i::RngIntElt :
    GenericFibreKinds := [], Minimalities := []) -> SeqEnum
{One model key per selected orbit over prefix i of FibrationPrefixes(R).
The parameter GenericFibreKinds (default []) selects "smooth", "regular cusp",
and "nonregular" when empty, or a subset of these case-sensitive values.
The parameter Minimalities (default []) selects "everywhere minimal" and
"nonminimal somewhere" when empty, or a subset of these case-sensitive values.
See also aInvariantsFromKey, WriteFibrationRepresentativeKeys}
    require valid_ring(R): "R must be a binary polynomial ring over GF(2) or GF(3)";
    P := prefixes(R);
    message := selection_error([i],GenericFibreKinds,Minimalities,#P);
    require message eq "": message;
    indices,kinds,minimalities := selections([i],GenericFibreKinds,Minimalities,#P);
    keys,cells := selected_representative_keys(R,P[i],kinds,minimalities);
    return keys;
end intrinsic;
