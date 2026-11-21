import "_archive_invariant.m": invariant_context, orbit_invariant;
import "_support.m": shell_quote;

// Canonical prefix labels allow the package to reorder or change representatives.
procedure test_archive(q, labels, expected)
    R<s, t> := PolynomialRing(GF(q), 2);
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
        if q eq 2 or label in labels then
            contexts[label] := ctx;
            indices[label] := i;
        end if;
    end for;
    if q eq 2 then labels := seen; assert #labels eq 88; end if;
    assert Keys(indices) eq labels;
    ordered := Sort(Setseq(labels));
    total, cells, rows := NumberOfFibrations(R : PrefixIndices := [indices[label] : label in ordered]);
    assert total eq expected[1][2];
    assert cells[<"smooth", "everywhere minimal">] eq expected[2][2];
    assert cells[<"smooth", "nonminimal somewhere">] eq expected[3][2];
    assert &+[cells[<kind, minimality>]
        : kind in ["regular cusp", "nonregular"],
          minimality in ["everywhere minimal", "nonminimal somewhere"]] eq expected[4][2];
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
            result := Pipe(Sprintf("python3 tests/support/_archive_hash.py prefix %o %o %o",
                q, label, shell_quote(output)), "");
            assert #Split(result, "\n") eq 4;
            combined cat:= result;
        end for;
        result := Pipe(Sprintf("python3 tests/support/_archive_hash.py combined %o", q), combined);
        assert result eq &cat[Sprintf("%o\t%o\t%o\n", row[1], row[2], row[3])
            : row in expected];
    catch e
        Remove(output);
        error e;
    end try;
    Remove(output);
end procedure;
