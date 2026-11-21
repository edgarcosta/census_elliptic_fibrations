import "tests/support/_support.m": compare_native;

procedure write_tables_reference(path)
    empty := Open(path, "w");
    delete empty;
    for p in [2, 3], d in [1..12] do
        fprintf path, "%o %o\n", p, d;
        fprintf path, "%o\n", Join([IntegerToString(Integers()!c)
            : c in Eltseq(ConwayPolynomial(p, d))], " ");
    end for;
end procedure;

compare_native(write_tables_reference, "tests/compare/tables.py");
