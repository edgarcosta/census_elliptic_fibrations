import "_archive_invariant.m": invariant_context, orbit_invariant;
import "_support.m": shell_quote;

// Canonical prefix labels allow the package to reorder or change representatives.
procedure test_archive(q)
    R<s, t> := PolynomialRing(GF(q), 2);
    expected := AssociativeArray();
    lines := Split(Read("tests/data/oldcode-archive/summaries.tsv"), "\n");
    for line in lines[2..#lines] do
        if line eq "CELLS" then break; end if;
        row := Split(line, "\t");
        if StringToInteger(row[1]) eq q then
            expected[<StringToInteger(row[3]), row[4]>] := StringToInteger(row[5]);
        end if;
    end for;
    labels := {key[1] : key in Keys(expected)};
    assert #labels eq (q eq 2 select 88 else 2);
    prefixes := FibrationPrefixes(R);
    contexts := AssociativeArray();
    indices := AssociativeArray();
    seen := {};
    for i in [1..#prefixes] do
        P := prefixes[i];
        W := q eq 2 select [R | P[1], 0, P[2], 0, 0]
            else [R | 0, P[1], 0, P[2], 0];
        ctx := invariant_context(W);
        label := ctx`Prefix;
        assert label notin seen;
        Include(~seen, label);
        if label in labels then
            contexts[label] := ctx;
            indices[label] := i;
        end if;
    end for;
    assert Keys(indices) eq labels;
    ordered := Sort(Setseq(labels));
    total, cells, rows := NumberOfFibrations(R : PrefixIndices := [indices[label] : label in ordered]);
    assert total eq &+[expected[<label, "all">] : label in ordered];
    for j in [1..#rows] do
        label := ordered[j];
        counts := rows[j][2];
        assert counts[<"smooth", "everywhere minimal">] eq expected[<label, "smooth_minimal">];
        assert counts[<"smooth", "nonminimal somewhere">] eq expected[<label, "smooth_nonminimal">];
        assert &+[counts[<kind, minimality>]
            : kind in ["regular cusp", "nonregular"],
              minimality in ["everywhere minimal", "nonminimal somewhere"]] eq expected[<label, "delta_zero">];
    end for;
    output := TemporaryFile();
    combined := "";
    try
        for label in ordered do
            O := Open(output, "wb");
            for cell in [1..3] do
                if cell eq 3 then
                    keys := FibrationRepresentativeKeys(R, indices[label] :
                        GenericFibreKinds := ["regular cusp", "nonregular"]);
                else
                    keys := FibrationRepresentativeKeys(R, indices[label] :
                        GenericFibreKinds := ["smooth"], Minimalities := [cell eq 1
                            select "everywhere minimal" else "nonminimal somewhere"]);
                end if;
                buffer := [Integers() | ];
                for key in keys do
                    W := aInvariantsFromKey(R, key);
                    invariant := orbit_invariant(W, contexts[label]);
                    buffer cat:= [(invariant div 256^b) mod 256 : b in [0..7]];
                    buffer cat:= [cell eq 3 select 0 else 1, cell eq 1 select 1 else 0];
                    if #buffer eq 10*4096 then
                        WriteBytes(O, buffer);
                        buffer := [Integers() | ];
                    end if;
                end for;
                if #buffer gt 0 then WriteBytes(O, buffer); end if;
            end for;
            Flush(O);
            delete O;
            result := Pipe(Sprintf("python3 tests/_archive_hash.py prefix %o %o %o",
                q, label, shell_quote(output)), "");
            assert #Split(result, "\n") eq 4;
            combined cat:= result;
        end for;
        result := Pipe(Sprintf("python3 tests/_archive_hash.py combined %o", q), combined);
        assert result eq "ARCHIVE_OK\n";
    catch e
        Remove(output);
        error e;
    end try;
    Remove(output);
end procedure;
