// Pins all six census cells against the reviewed literal census totals.
procedure test_six_cell_counts_f2()
    R<s,t> := PolynomialRing(GF(2),2);
    assert #FibrationPrefixes(R) eq 88;
    n,cells,rows := NumberOfFibrations(R);
    assert n eq 1402468 and #rows eq 88;
    expected := [1397598, 4134, 680, 38, 7, 11];
    j := 0;
    for kind in ["smooth","regular cusp","nonregular"] do
        for minimality in ["everywhere minimal","nonminimal somewhere"] do
            j +:= 1;
            cell := <kind,minimality>;
            assert cells[cell] eq expected[j];
            assert &+[ row[2][cell] : row in rows ] eq expected[j];
        end for;
    end for;
    for kind in ["smooth","regular cusp","nonregular"] do
        for minimality in ["everywhere minimal","nonminimal somewhere"] do
            selected,table,part := NumberOfFibrations(R : PrefixIndices := [2,1],
                GenericFibreKinds := [kind], Minimalities := [minimality]);
            assert [ row[1] : row in part ] eq [2,1];
            assert #Keys(table) eq 1;
            assert selected eq rows[1][2][<kind,minimality>]+rows[2][2][<kind,minimality>];
        end for;
    end for;
    for kind in ["smooth","regular cusp","nonregular"] do
        selected,table,part := NumberOfFibrations(R : PrefixIndices := [1], GenericFibreKinds := [kind]);
        assert #Keys(table) eq 2;
        assert selected eq &+[ table[cell] : cell in Keys(table) ];
    end for;
    for minimality in ["everywhere minimal","nonminimal somewhere"] do
        selected,table,part := NumberOfFibrations(R : PrefixIndices := [1], Minimalities := [minimality]);
        assert #Keys(table) eq 3;
        assert selected eq &+[ table[cell] : cell in Keys(table) ];
    end for;
end procedure;

test_six_cell_counts_f2();
