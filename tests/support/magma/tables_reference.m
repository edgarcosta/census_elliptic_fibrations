SetQuitOnError(true);
AttachSpec("spec");

empty := Open(referencepath, "w");
delete empty;
for p in [2, 3], d in [1..12] do
    table := FibrationTraceTable(p, d);
    fprintf referencepath, "%o %o\n", p, d;
    fprintf referencepath, "%o\n",
        Join([IntegerToString(Integers()!c) : c in Eltseq(table`Polynomial)], " ");
    ordinary := table`Ordinary;
    if d le 4 then
        K := table`Field;
        ordinary := [TraceOfFrobenius(EllipticCurve(p eq 2
            select [K | 1, 0, 0, 0, table`Generator^e]
            else [K | 0, 1, 0, 0, table`Generator^e])) : e in [0..p^d-2]];
    end if;
    fprintf referencepath, "%o\n", Join([IntegerToString(t) : t in ordinary], " ");
    assert table`NodalTrace eq 1;
end for;
quit;
