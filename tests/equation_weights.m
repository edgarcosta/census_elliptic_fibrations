import "tests/_weierstrass_change.m": homogeneous_forms;

// Checks each raw prefix orbit and translation fibre independently of the engine.
procedure test_equation_weights(p)
    k := GF(p);
    R<s, t> := PolynomialRing(k, 2);
    S<x, y> := PolynomialRing(R, 2);
    shifts := homogeneous_forms(R, p eq 2 select 2 else 4);
    lower_shifts := homogeneous_forms(R, p eq 2 select 6 else 4);
    base := [Matrix(g) : g in GL(2, k)];
    prefixes := FibrationPrefixes(R);
    assert #prefixes eq (p eq 2 select 88 else 1860);
    total, _, rows := NumberOfEquations(R);
    assert [row[1] : row in rows] eq [1..#prefixes];
    dimension := p eq 2 select 22 else 13;
    seen := {};
    raw_total := 0;
    for i in [1..#prefixes] do
        a, b := Explode(prefixes[i]);
        equation := p eq 2 select y^2+a*x*y+b*y-x^3
            else y^2-x^3-a*x^2-b*x;
        prefix_image := {R!0};
        images := {};
        for v in lower_shifts do
            changed := p eq 2 select Evaluate(equation, [x, y+v])
                else Evaluate(equation, [x+v, y]);
            difference := changed-equation;
            if p eq 3 then
                delta := -MonomialCoefficient(difference, x);
                Include(~prefix_image, delta);
                if delta ne 0 then continue; end if;
            end if;
            Include(~images, difference);
        end for;
        orbit := {};
        stabilizer := 0;
        for g in base do
            coordinates := [g[1, 1]*s+g[1, 2]*t, g[2, 1]*s+g[2, 2]*t];
            ag := Evaluate(a, coordinates);
            bg := Evaluate(b, coordinates);
            if p eq 2 then
                equation := y^2+ag*x*y+bg*y-x^3;
                for u in shifts do
                    changed := Evaluate(equation, [x+u^2+ag*u, y+u*x]);
                    assert MonomialCoefficient(changed, x^2) eq 0;
                    pair := <MonomialCoefficient(changed, x*y),
                        MonomialCoefficient(changed, y)>;
                    Include(~orbit, pair);
                    if pair eq prefixes[i] then stabilizer +:= 1; end if;
                end for;
            else
                // A remainder labels a coset of the quartic x-translations.
                remainder := ag eq 0 select bg else NormalForm(bg, [ag]);
                Include(~orbit, <ag, remainder>);
                if ag eq a and b-bg in prefix_image then
                    stabilizer +:= #shifts div #prefix_image;
                end if;
            end if;
        end for;
        raw_size := #orbit*#prefix_image;
        assert raw_size*stabilizer eq #base*#shifts;
        error if #(seen meet orbit) ne 0,
            Sprintf("GF(%o), prefix %o: repeated raw prefix orbit", p, i);
        seen join:= orbit;
        raw_total +:= raw_size;
        weight := raw_size*#images;
        cells := rows[i][2];
        error if &+[cells[cell] : cell in Keys(cells)] ne raw_size*p^dimension,
            Sprintf("GF(%o), prefix %o: wrong equation-row size", p, i);
        error if &or[cells[cell] mod weight ne 0 : cell in Keys(cells)],
            Sprintf("GF(%o), prefix %o: nonintegral translation-fibre count", p, i);
    end for;
    assert raw_total eq (p eq 2 select 2^10 else 3^14);
    assert total eq (p eq 2 select 2^32 else 3^27);
end procedure;

test_equation_weights(2);
test_equation_weights(3);
