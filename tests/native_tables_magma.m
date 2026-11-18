import "tests/support/_support.m": expect_error;

procedure test_tables()
    directory := "tests/fixtures/trace_tables/";
    hashes := AssociativeArray();
    for line in Split(Read(directory cat "tables.sha256"), "\n") do
        parts := Split(line, " ");
        hashes[parts[2]] := parts[1];
    end for;
    for p in [2, 3], d in [1..12] do
        filename := Sprintf("f%o-d%o.txt", p, d);
        lines := Split(Read(directory cat filename), "\n");
        assert lines[1] eq Sprintf("%o %o", p, d);
        polynomial := [StringToInteger(n) : n in Split(lines[2], " ")];
        legacy := [StringToInteger(n) : n in Split(lines[3], " ")];
        table := FibrationTraceTable(p, d);
        K := table`Field;
        alpha := table`Generator;
        assert #K eq p^d;
        assert [Integers()!c : c in Eltseq(table`Polynomial)] eq polynomial;
        assert Evaluate(table`Polynomial, alpha) eq 0 and Order(alpha) eq p^d-1;
        assert #table`Ordinary eq p^d-1 and table`Ordinary eq legacy;
        assert table`NodalTrace eq 1;
        assert &and[Abs(t) le p^d : t in table`Ordinary];
        if d le 4 then
            a := K!1;
            for e in [0..p^d-2] do
                E := EllipticCurve(p eq 2 select [K | 1, 0, 0, 0, a]
                    else [K | 0, 1, 0, 0, a]);
                assert TraceOfFrobenius(E) eq table`Ordinary[e+1];
                a *:= alpha;
            end for;
        end if;
        canonical := Sprintf("%o %o\n", p, d)
            cat Join([IntegerToString(n) : n in polynomial], " ") cat "\n"
            cat Join([IntegerToString(n) : n in table`Ordinary], " ") cat "\n";
        digest := Split(Pipe("sha256sum", canonical), " ")[1];
        assert digest eq hashes[filename];
    end for;
    for pair in [<2, 0>, <3, 13>, <5, 1>] do
        expect_error(procedure() t := FibrationTraceTable(pair[1], pair[2]); end procedure,
            "p must be 2 or 3 and degree must be between 1 and 12");
    end for;
end procedure;

test_tables();
