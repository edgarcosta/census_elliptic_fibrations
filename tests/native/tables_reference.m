SetQuitOnError(true);
AttachSpec("spec");

for p in [2, 3], d in [1..12] do
    table := FibrationTraceTable(p, d);
    fprintf referencepath, "%o %o\n", p, d;
    fprintf referencepath, "%o\n",
        Join([IntegerToString(Integers()!c) : c in Eltseq(table`Polynomial)], " ");
    fprintf referencepath, "%o\n", Join([IntegerToString(t) : t in table`Ordinary], " ");
    assert table`NodalTrace eq 1;
end for;
quit;
