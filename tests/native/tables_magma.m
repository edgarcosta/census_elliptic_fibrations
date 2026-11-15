SetQuitOnError(true);
AttachSpec("spec");
import "tests/_support.m": expect_error;

// Counting both coordinates catches trace signs without a character or transform oracle.
function direct_trace(K, a)
    p := Characteristic(K);
    affine := 0;
    for x, y in K do
        if p eq 2 then
            affine +:= y^2+x*y eq x^3+a select 1 else 0;
        else
            affine +:= y^2 eq x^3+x^2+a select 1 else 0;
        end if;
    end for;
    return #K-affine;
end function;

procedure test_tables()
    directory := "tests/native/data/trace_tables/";
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
            assert direct_trace(K, K!0) eq table`NodalTrace;
            a := K!1;
            for e in [0..p^d-2] do
                assert direct_trace(K, a) eq table`Ordinary[e+1];
                a *:= alpha;
            end for;
        end if;
        canonical := Sprintf("%o %o\n", p, d)
            cat Join([IntegerToString(n) : n in polynomial], " ") cat "\n"
            cat Join([IntegerToString(n) : n in table`Ordinary], " ") cat "\n";
        digest := Split(Pipe("sha256sum", canonical), " ")[1];
        assert digest eq hashes[filename];
        printf "%o  %o\n", digest, filename;
    end for;
    for pair in [<2, 0>, <3, 13>, <5, 1>] do
        expect_error(procedure() t := FibrationTraceTable(pair[1], pair[2]); end procedure,
            "p must be 2 or 3 and degree must be between 1 and 12");
    end for;
end procedure;

test_tables();
Write(resultpath, "TABLES_MAGMA_OK");
quit;
