import "Tests/_support.m": expect_error;

// Pins all six equation cells against local-density and singular-locus formulas.
procedure test_equation_closed_forms(p)
    R<s,t> := PolynomialRing(GF(p),2);
    if p eq 2 then
        all_equations := 2^(3+7+9+13);
        cusps := 2^22-2^12;
        nonregular := 2^15;
        // Three rational places and one quadratic, minus three rational pairs.
        nonminimal := 3*2^22-2*2^12;
        cusp_nonminimal := 3*(2^16-2^11);
        nonregular_nonminimal := 3*2^13-2*2^11;
    else
        all_equations := 3^(5+9+13);
        cusps := 3^13-3^5;
        nonregular := 3^10;
        // Four rational places and three quadratics, minus six rational pairs.
        nonminimal := 4*3^17-3^8;
        cusp_nonminimal := 4*(3^9-3^5);
        nonregular_nonminimal := (1+4*(3^3-1)+(3-6)*(3-1))*3^5;
    end if;
    smooth_nonminimal := nonminimal-cusp_nonminimal-nonregular_nonminimal;
    expected := [ all_equations-cusps-nonregular-smooth_nonminimal,
        smooth_nonminimal, cusps-cusp_nonminimal, cusp_nonminimal,
        nonregular-nonregular_nonminimal, nonregular_nonminimal ];
    n,cells,rows := NumberOfEquations(R);
    assert n eq all_equations;
    assert #Keys(cells) eq 6;
    assert [ row[1] : row in rows ] eq [1..#FibrationPrefixes(R)];
    j := 0;
    for kind in ["smooth","regular cusp","nonregular"] do
        for minimality in ["everywhere minimal","nonminimal somewhere"] do
            j +:= 1;
            cell := <kind,minimality>;
            error if cells[cell] ne expected[j], Sprintf(
                "GF(%o), %o: engine %o, closed form %o",p,cell,cells[cell],expected[j]);
            assert &+[ row[2][cell] : row in rows ] eq expected[j];
        end for;
    end for;
    selected,table,part := NumberOfEquations(R : PrefixIndices := [2,1],
        GenericFibreKinds := ["nonregular","smooth"],
        Minimalities := ["nonminimal somewhere"]);
    assert [ row[1] : row in part ] eq [2,1];
    assert Keys(table) eq { <kind,"nonminimal somewhere"> : kind in ["smooth","nonregular"] };
    assert selected eq &+[ table[cell] : cell in Keys(table) ];
    for cell in Keys(table) do
        assert table[cell] eq rows[1][2][cell]+rows[2][2][cell];
        assert part[1][2][cell] eq rows[2][2][cell];
        assert part[2][2][cell] eq rows[1][2][cell];
    end for;
end procedure;

// Pins the same rejected selections and rings as the orbit-count interface.
procedure test_equation_count_errors()
    R := PolynomialRing(GF(2),2);
    expect_error(procedure() n := NumberOfEquations(R : PrefixIndices := {1}); end procedure,
        "Selections must be sequences");
    expect_error(procedure() n := NumberOfEquations(R : PrefixIndices := [0]); end procedure,
        "valid prefix indices");
    expect_error(procedure() n := NumberOfEquations(R : PrefixIndices := [1,1]); end procedure,
        "PrefixIndices must be distinct");
    expect_error(procedure() n := NumberOfEquations(R : GenericFibreKinds := ["cusp"]); end procedure,
        "GenericFibreKinds must contain distinct values");
    expect_error(procedure() n := NumberOfEquations(R : Minimalities := ["minimal"]); end procedure,
        "Minimalities must contain distinct values");
    S := PolynomialRing(GF(2^2),2);
    expect_error(procedure() n := NumberOfEquations(S); end procedure,"binary polynomial ring");
    S := PolynomialRing(GF(2),3);
    expect_error(procedure() n := NumberOfEquations(S); end procedure,"binary polynomial ring");
end procedure;

test_equation_closed_forms(2);
test_equation_closed_forms(3);
test_equation_count_errors();
