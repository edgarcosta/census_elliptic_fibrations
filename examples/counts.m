// Verifies the public census totals over both supported fields.
procedure example_counts()
    for p in [2,3] do
        R := PolynomialRing(GF(p),2);
        n,cells,rows := NumberOfFibrations(R);
        assert n eq (p eq 2 select 1402468 else 1307814237);
        assert #rows gt 0;
    end for;
end procedure;

example_counts();
